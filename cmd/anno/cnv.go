package anno

import (
	"log"
	"open-anno/anno"
	"open-anno/pkg"
	"os"
	"path"
	"strings"
	"sync"

	"github.com/brentp/bix"
	"github.com/brentp/vcfgo"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type AnnoCnvParam struct {
	Input              string   `validate:"required,pathexists"`
	GenePred           string   `validate:"required,pathexists"`
	GenePredIndex      string   `validate:"required,pathexists"`
	Gene               string   `validate:"required,pathexists"`
	Output             string   `validate:"required"`
	RegionBaseds       []string `validate:"pathsexists"`
	RegionBasedIndexes []string `validate:"pathsexists"`
	Overlap            float64  `validate:"required"`
	Concurrency        int      `validate:"required"`
}

func (this *AnnoCnvParam) Valid() error {
	this.GenePredIndex = this.GenePred + ".tbi"
	for _, db := range this.RegionBaseds {
		this.RegionBasedIndexes = append(this.RegionBasedIndexes, db+".tbi")
	}
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	validate.RegisterValidation("pathsexists", pkg.CheckPathsExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return os.MkdirAll(this.Outdir(), 0666)
}

func (this AnnoCnvParam) Outdir() string {
	return path.Dir(this.Output)
}

func (this *AnnoCnvParam) RunAnno(cnvs []*pkg.CNV, gpeTbx *bix.Bix, rbTbxs []*bix.Bix, dbnames []string) (map[string]map[string]any, error) {
	cnvChan := make(chan *pkg.CNV, len(cnvs))
	for _, snv := range cnvs {
		cnvChan <- snv
	}
	close(cnvChan)
	var wg sync.WaitGroup
	resChan := make(chan anno.AnnoInfo, len(cnvs))
	for i := 0; i <= this.Concurrency; i++ {
		wg.Add(1)
		go anno.AnnoCnvWorker(cnvChan, gpeTbx, rbTbxs, dbnames, this.Overlap, resChan, &wg)
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

func (this AnnoCnvParam) Run() error {
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
	// 打开GenePred
	gpeTbx, err := bix.New(this.GenePred)
	if err != nil {
		return err
	}
	defer gpeTbx.Close()
	vcfHeader.Infos["DETAIL"] = &vcfgo.Info{
		Id:          "DETAIL",
		Description: "Gene detail, FORMAT=Gene:Transcript:Exon:NA_CHANGE:AA_CHANGE",
		Number:      ".",
		Type:        "String",
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
	// 读取变异
	cnvs := make([]*pkg.CNV, 0)
	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
		if len(variant.Chrom()) > 5 {
			continue
		}
		cnvs = append(cnvs, &pkg.CNV{Variant: *variant})
	}
	annoResult, err := this.RunAnno(cnvs, gpeTbx, rbTbxs, dbnames)
	if err != nil {
		return err
	}
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	vcfWriter, err := vcfgo.NewWriter(writer, vcfHeader)
	for _, cnv := range cnvs {
		for key, val := range annoResult[cnv.PK()] {
			if val != "" && val != "." {
				err = cnv.Info().Set(key, val)
				if err != nil {
					return err
				}
			}
		}
		vcfWriter.WriteVariant(&cnv.Variant)
	}
	return nil
}

func NewAnnoCnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "cnv",
		Short: "Annotate CNV",
		Run: func(cmd *cobra.Command, args []string) {
			var param AnnoCnvParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.Gene, _ = cmd.Flags().GetString("gene")
			param.Output, _ = cmd.Flags().GetString("output")
			param.RegionBaseds, _ = cmd.Flags().GetStringArray("regionbaseds")
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
	cmd.Flags().StringP("gene", "g", "", "Input Gene Symbol To ID File")
	cmd.Flags().StringP("gbname", "n", "", "Parameter Database Name")
	cmd.Flags().StringP("output", "o", "", "AnnoOutput File")
	cmd.Flags().StringArrayP("regionbaseds", "r", []string{}, "Input RegionBased Database File")
	cmd.Flags().Float64P("overlap", "l", 0.7, "Parameter Database Name")
	cmd.Flags().IntP("concurrency", "c", 10000, "Parameter Concurrency Numbers")
	return cmd
}
