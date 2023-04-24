package anno

import (
	"log"
	"open-anno/anno"
	"open-anno/anno/db"
	"open-anno/pkg"
	"os"
	"path"
	"sort"
	"strings"
	"sync"

	"github.com/brentp/bix"
	"github.com/brentp/faidx"
	"github.com/brentp/vcfgo"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
	"github.com/syndtr/goleveldb/leveldb"
	"github.com/syndtr/goleveldb/leveldb/opt"
)

type AnnoSnvParam struct {
	Input               string `validate:"required,pathexists"`
	GenePred            string `validate:"required,pathexists"`
	GenePredIndex       string `validate:"required,pathexists"`
	Genome              string `validate:"required,pathexists"`
	GenomeIndex         string `validate:"required,pathexists"`
	Gene                string `validate:"required,pathexists"`
	Output              string `validate:"required"`
	AAshort             bool
	Exon                bool
	FilterBaseds        []string `validate:"pathsexists"`
	FilterBasedIndexes  []string `validate:"pathsexists"`
	RegionBaseds        []string `validate:"pathsexists"`
	RegionBasedIndexes  []string `validate:"pathsexists"`
	FilterBasedLevelDBs []string `validate:"pathsexists"`
	Overlap             float64  `validate:"required"`
	Concurrency         int      `validate:"required"`
}

func (this *AnnoSnvParam) Valid() error {
	this.GenePredIndex = this.GenePred + ".tbi"
	this.GenomeIndex = this.Genome + ".fai"
	for _, db := range this.FilterBaseds {
		this.FilterBasedIndexes = append(this.FilterBasedIndexes, db+".tbi")
	}
	for _, db := range this.RegionBaseds {
		this.RegionBasedIndexes = append(this.RegionBasedIndexes, db+".tbi")
	}
	pkg.IS_EXON_REGION = this.Exon
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	validate.RegisterValidation("pathsexists", pkg.CheckPathsExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return os.MkdirAll(this.Outdir(), 0666)
}

func (this AnnoSnvParam) Outdir() string {
	return path.Dir(this.Output)
}

func (this AnnoSnvParam) RunAnno(snvs []*pkg.SNV, gpeTbx *bix.Bix, fbTbxs []*bix.Bix, rbTbxs []*bix.Bix, dbnames []string, fbDBs []*leveldb.DB, genome *faidx.Faidx) (map[string]map[string]any, error) {
	snvChan := make(chan *pkg.SNV, len(snvs))
	for _, snv := range snvs {
		snvChan <- snv
	}
	close(snvChan)
	var wg sync.WaitGroup
	resChan := make(chan anno.AnnoInfo, len(snvs))
	for i := 0; i <= this.Concurrency; i++ {
		wg.Add(1)
		go anno.AnnoSnvWorker(snvChan, gpeTbx, fbTbxs, rbTbxs, dbnames, fbDBs, genome, this.Overlap, resChan, &wg)
	}
	go func() {
		wg.Wait()
		close(resChan)
	}()
	results := make(map[string]map[string]any)
	for res := range resChan {
		if res.Error != nil {
			return results, res.Error
		}
		results[res.PK] = res.Data
	}
	return results, nil
}

func (this AnnoSnvParam) Run() error {
	// 读取GeneID信息
	log.Printf("Read Gene: %s ...", this.Gene)
	err := pkg.InitGeneSymbolToID(this.Gene)
	if err != nil {
		return err
	}
	// 打开变异输入文件
	log.Printf("Read AnnoInput: %s ...", this.Input)
	reader, err := pkg.NewIOReader(this.Input)
	if err != nil {
		return err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return err
	}
	defer vcfReader.Close()
	vcfHeader := vcfReader.Header
	snvs := make([]*pkg.SNV, 0)
	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
		if len(variant.Chrom()) > 5 {
			continue
		}
		snvs = append(snvs, &pkg.SNV{Variant: *variant})
	}
	// 打开GenePred
	log.Printf("Open TABIX Handle ...")
	gpeTbx, err := bix.New(this.GenePred)
	if err != nil {
		return err
	}
	defer gpeTbx.Close()
	vcfHeader.Infos["GENE"] = &vcfgo.Info{
		Id:          "GENE",
		Description: "Gene Symbol",
		Number:      ".",
		Type:        "String",
	}
	vcfHeader.Infos["GENE_ID"] = &vcfgo.Info{
		Id:          "GENE_ID",
		Description: "Gene Entrez ID",
		Number:      ".",
		Type:        "String",
	}
	vcfHeader.Infos["REGION"] = &vcfgo.Info{
		Id:          "REGION",
		Description: "Region in gene, eg: exonic, intronic, UTR3, UTR5",
		Number:      ".",
		Type:        "String",
	}
	vcfHeader.Infos["EVENT"] = &vcfgo.Info{
		Id:          "EVENT",
		Description: "Variant Event, eg: missense, nonsense, splicing",
		Number:      ".",
		Type:        "String",
	}
	vcfHeader.Infos["DETAIL"] = &vcfgo.Info{
		Id:          "DETAIL",
		Description: "Gene detail, FORMAT=Gene:Transcript:Exon:NA_CHANGE:AA_CHANGE",
		Number:      ".",
		Type:        "String",
	}
	// 打开FilterBaseds
	fbTbxs := make([]*bix.Bix, len(this.FilterBaseds))
	for i, fbFile := range this.FilterBaseds {
		fbTbxs[i], err = bix.New(fbFile)
		if err != nil {
			return err
		}
		for key, info := range fbTbxs[i].VReader.Header.Infos {
			vcfHeader.Infos[key] = info
		}
		defer fbTbxs[i].Close()
	}
	// 打开RegionBaseds
	rbTbxs := make([]*bix.Bix, len(this.RegionBaseds))
	dbnames := make([]string, len(this.RegionBaseds))
	for i, rbFile := range this.RegionBaseds {
		rbTbxs[i], err = bix.New(rbFile)
		if err != nil {
			return err
		}
		dbnames[i] = strings.Split(path.Base(rbFile), ".")[0]
		vcfHeader.Infos[dbnames[i]] = &vcfgo.Info{
			Id:          dbnames[i],
			Description: dbnames[i],
			Number:      ".",
			Type:        "String",
		}
		defer rbTbxs[i].Close()
	}
	// 打开FilterBasedLevelDBs
	log.Printf("Open LevelDB ...")
	fbDBs := make([]*leveldb.DB, len(this.FilterBasedLevelDBs))
	for i, fbDir := range this.FilterBasedLevelDBs {
		fbDBs[i], err = leveldb.OpenFile(fbDir, &opt.Options{ReadOnly: true})
		if err != nil {
			return err
		}
		headerInfo, err := db.GetHeaderLevelDB(fbDBs[i])
		if err != nil {
			return err
		}
		for _, info := range headerInfo {
			vcfHeader.Infos[info.Id] = info
		}
		defer fbDBs[i].Close()
	}
	// 打开Genome
	genome, err := faidx.New(this.Genome)
	if err != nil {
		return err
	}
	annoResult, err := this.RunAnno(snvs, gpeTbx, fbTbxs, rbTbxs, dbnames, fbDBs, genome)
	if err != nil {
		return err
	}
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	vcfWriter, err := vcfgo.NewWriter(writer, vcfHeader)
	whiteList := []string{"GENE", "GENE_ID", "EVENT", "REGION", "DETAIL"}
	for _, snv := range snvs {
		for key, val := range annoResult[snv.PK()] {
			idx := sort.SearchStrings(whiteList, key)
			if (idx < len(whiteList) && whiteList[idx] == key) || (val != "" && val != ".") {
				err = snv.Info().Set(key, val)
				if err != nil {
					return err
				}
			}
		}
		vcfWriter.WriteVariant(&snv.Variant)
	}
	return nil
}

func NewAnnoSnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "snv",
		Short: "Annotate  SNV",
		Run: func(cmd *cobra.Command, args []string) {
			var param AnnoSnvParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.Genome, _ = cmd.Flags().GetString("genome")
			param.Gene, _ = cmd.Flags().GetString("gene")
			param.Output, _ = cmd.Flags().GetString("output")
			param.AAshort, _ = cmd.Flags().GetBool("aashort")
			param.Exon, _ = cmd.Flags().GetBool("exon")
			param.FilterBaseds, _ = cmd.Flags().GetStringArray("filterbaseds")
			param.RegionBaseds, _ = cmd.Flags().GetStringArray("regionbaseds")
			param.FilterBasedLevelDBs, _ = cmd.Flags().GetStringArray("leveldbs")
			param.Overlap, _ = cmd.Flags().GetFloat64("overlap")
			param.Concurrency, _ = cmd.Flags().GetInt("concurrency")
			err := param.Valid()
			if err != nil {
				cmd.Help()
				log.Fatal(err)
			}
			err = param.Run()
			if err != nil {
				log.Fatal(err)
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "AnnoInput File")
	cmd.Flags().StringP("genepred", "d", "", "Input GenePred File")
	cmd.Flags().StringP("genome", "G", "", "Input Genome Fasta File")
	cmd.Flags().StringP("gene", "g", "", "Input Gene Symbol To ID File")
	cmd.Flags().StringP("output", "o", "", "AnnoOutput File")
	cmd.Flags().BoolP("aashort", "a", false, "Parameter Is AA Short")
	cmd.Flags().BoolP("exon", "e", false, "Parameter Is Exon")
	cmd.Flags().StringArrayP("filterbaseds", "f", []string{}, "Input FilterBased Database File")
	cmd.Flags().StringArrayP("regionbaseds", "r", []string{}, "Input RegionBased Database File")
	cmd.Flags().StringArrayP("leveldbs", "F", []string{}, "Input FilterBased LevelDB Directory")
	cmd.Flags().Float64P("overlap", "l", 0.7, "Parameter Database Name")
	cmd.Flags().IntP("concurrency", "c", 4, "Parameter Concurrency Numbers")
	return cmd
}
