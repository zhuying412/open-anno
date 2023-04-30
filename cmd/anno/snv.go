package anno

import (
	"fmt"
	"io/ioutil"
	"log"
	"open-anno/anno"
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
)

type AnnoSnvParam struct {
	Input              string `validate:"required,pathexists"`
	GenePred           string `validate:"required,pathexists"`
	GenePredIndex      string `validate:"required,pathexists"`
	Genome             string `validate:"required,pathexists"`
	GenomeIndex        string `validate:"required,pathexists"`
	Gene               string `validate:"required,pathexists"`
	Output             string `validate:"required"`
	AAshort            bool
	Exon               bool
	FilterBaseds       []string `validate:"pathsexists"`
	FilterBasedIndexes []string `validate:"pathsexists"`
	RegionBaseds       []string `validate:"pathsexists"`
	RegionBasedIndexes []string `validate:"pathsexists"`
	FilterBasedDirs    []string `validate:"pathsexists"`
	Overlap            float64  `validate:"required"`
	Concurrency        int      `validate:"required"`
	Chrom              string
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

func (this AnnoSnvParam) RunAnno(snvs []*pkg.SNV, gpeTbx *bix.Bix, fbTbxs []*bix.Bix, rbTbxs []*bix.Bix, dbnames []string, genome *faidx.Faidx) (map[string]map[string]any, error) {
	snvChan := make(chan *pkg.SNV, len(snvs))
	for _, snv := range snvs {
		snvChan <- snv
	}
	close(snvChan)
	var wg sync.WaitGroup
	resChan := make(chan anno.AnnoInfo, len(snvs))
	for i := 0; i <= this.Concurrency; i++ {
		wg.Add(1)
		go anno.AnnoSnvWorker(snvChan, gpeTbx, fbTbxs, rbTbxs, dbnames, genome, this.Overlap, resChan, &wg)
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
func (this AnnoSnvParam) GetHeaderInfos() (map[string]*vcfgo.Info, []string, error) {
	infos := map[string]*vcfgo.Info{
		"GENE":    {Id: "GENE", Description: "Gene Symbol", Number: ".", Type: "String"},
		"GENE_ID": {Id: "GENE_ID", Description: "Gene Entrez ID", Number: ".", Type: "String"},
		"REGION":  {Id: "REGION", Description: "Region in gene, eg: exonic, intronic, UTR3, UTR5", Number: ".", Type: "String"},
		"EVENT":   {Id: "EVENT", Description: "Variant Event, eg: missense, nonsense, splicing", Number: ".", Type: "String"},
		"DETAIL":  {Id: "DETAIL", Description: "Gene detail, FORMAT=Gene:Transcript:Exon:NA_CHANGE:AA_CHANGE", Number: ".", Type: "String"},
	}
	for _, fbFile := range this.FilterBaseds {
		fbTbx, err := bix.New(fbFile)
		if err != nil {
			return infos, []string{}, err
		}
		defer fbTbx.Close()
		for key, info := range fbTbx.VReader.Header.Infos {
			infos[key] = info
		}
	}
	// 读取目录下的文件和子目录
	for _, fbDir := range this.FilterBasedDirs {
		fbFiles, err := ioutil.ReadDir(fbDir)
		if err != nil {
			return infos, []string{}, err
		}
		for _, fbFile := range fbFiles {
			if strings.HasSuffix(fbFile.Name(), ".vcf.gz") {
				fbTbx, err := bix.New(path.Join(fbDir, fbFile.Name()))
				if err != nil {
					return infos, []string{}, err
				}
				defer fbTbx.Close()
				for key, info := range fbTbx.VReader.Header.Infos {
					infos[key] = info
				}
				break
			}
		}
	}
	dbnames := make([]string, len(this.RegionBaseds))
	for i, rbFile := range this.RegionBaseds {
		dbname := strings.Split(path.Base(rbFile), ".")[0]
		dbnames[i] = dbname
		infos[dbname] = &vcfgo.Info{Id: dbname, Description: dbname, Number: ".", Type: "String"}
	}
	return infos, dbnames, nil
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
	keys := make([]string, 0)
	snvsMap := make(map[string][]*pkg.SNV)
	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
		snv := &pkg.SNV{Variant: *variant}
		chrom := variant.Chrom()
		if len(chrom) > 5 || (this.Chrom != "" && chrom != this.Chrom) {
			continue
		}
		key := fmt.Sprintf("%s.%d", chrom, snv.Pos/uint64(pkg.FilterBasedBucketSize))
		if snvs, ok := snvsMap[key]; !ok {
			snvsMap[key] = []*pkg.SNV{snv}
			keys = append(keys, key)
		} else {
			snvsMap[key] = append(snvs, snv)
		}
	}
	// VcfHeaderInfo
	infos, dbnames, err := this.GetHeaderInfos()
	if err != nil {
		return err
	}
	for id, info := range infos {
		vcfHeader.Infos[id] = info
	}
	// 打开Genome
	log.Printf("Open Genome Faidx ...")
	genome, err := faidx.New(this.Genome)
	if err != nil {
		return err
	}
	defer genome.Close()
	// 打开GenePred
	log.Printf("Open TABIX Handle ...")
	gpeTbx, err := bix.New(this.GenePred)
	if err != nil {
		return err
	}
	defer gpeTbx.Close()
	// 打开FilterBaseds
	fbTbxs := make([]*bix.Bix, len(this.FilterBaseds))
	for i, fbFile := range this.FilterBaseds {
		fbTbxs[i], err = bix.New(fbFile)
		if err != nil {
			return err
		}
		defer fbTbxs[i].Close()
	}
	// 打开RegionBaseds
	rbTbxs := make([]*bix.Bix, len(this.RegionBaseds))
	for i, rbFile := range this.RegionBaseds {
		rbTbxs[i], err = bix.New(rbFile)
		if err != nil {
			return err
		}
		defer rbTbxs[i].Close()
	}
	// 打开输出句柄
	log.Printf("Write to %s ...", this.Output)
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	vcfWriter, err := vcfgo.NewWriter(writer, vcfHeader)
	// 开始注释
	log.Println("Run Annotating ...")
	whiteList := []string{"GENE", "GENE_ID", "EVENT", "REGION", "DETAIL"}
	for _, key := range keys {
		log.Printf("Run Annotating %s ...", key)
		snvs, ok := snvsMap[key]
		if !ok {
			continue
		}
		fbTbxs2 := make([]*bix.Bix, 0)
		for _, fbDir := range this.FilterBasedDirs {
			fbFile := fmt.Sprintf("%s/%s.vcf.gz", fbDir, key)
			_, err := os.Stat(fbFile)
			if os.IsNotExist(err) {
				continue
			}
			fbTbx, err := bix.New(fbFile)
			if err != nil {
				return err
			}
			fbTbxs2 = append(fbTbxs2, fbTbx)
		}
		annoResult, err := this.RunAnno(snvs, gpeTbx, append(fbTbxs, fbTbxs2...), rbTbxs, dbnames, genome)
		for _, fbTbx := range fbTbxs2 {
			fbTbx.Close()
		}
		if err != nil {
			return err
		}
		for _, snv := range snvs {
			for id, val := range annoResult[snv.PK()] {
				idx := sort.SearchStrings(whiteList, id)
				if (idx < len(whiteList) && whiteList[idx] == id) || (val != "" && val != ".") {
					err = snv.Info().Set(id, val)
					if err != nil {
						return err
					}
				}
			}
			vcfWriter.WriteVariant(&snv.Variant)
		}
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
			param.FilterBasedDirs, _ = cmd.Flags().GetStringArray("filterbased_dirs")
			param.Overlap, _ = cmd.Flags().GetFloat64("overlap")
			param.Concurrency, _ = cmd.Flags().GetInt("concurrency")
			param.Chrom, _ = cmd.Flags().GetString("chrom")
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
	cmd.Flags().StringArrayP("filterbased_dirs", "F", []string{}, "Input FilterBased Directory")
	cmd.Flags().Float64P("overlap", "l", 0.7, "Parameter Database Name")
	cmd.Flags().IntP("concurrency", "c", 4, "Parameter Concurrency Numbers")
	cmd.Flags().StringP("chrom", "m", "", "Chromosome")
	return cmd
}
