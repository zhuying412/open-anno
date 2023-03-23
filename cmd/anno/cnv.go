package anno

import (
	"log"
	"open-anno/anno"
	"open-anno/anno/db"
	"open-anno/anno/gene"
	"open-anno/pkg"
	"os"
	"path"

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
	Goroutine          int      `validate:"required"`
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

func (this *AnnoCnvParam) RunAnnoGeneBase() (anno.AnnoResult, error) {
	// 读取gene输入文件
	log.Printf("Read Gene: %s ...", this.Gene)
	err := pkg.InitGeneSymbolToID(this.Gene)
	if err != nil {
		return anno.AnnoResult{}, err
	}
	return gene.AnnoCnvs(this.Input, this.GenePred, this.Goroutine)
}

func (this AnnoCnvParam) Run() error {
	// 读取变异输入文件
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
	// 构造 error channel
	var annoResults []anno.AnnoResult
	var annoResult anno.AnnoResult
	annoResult, err = this.RunAnnoGeneBase()
	if err != nil {
		return err
	}
	annoResults = append(annoResults, annoResult)
	// RegionBased 注释
	for _, database := range this.RegionBaseds {
		annoResult, err = db.AnnoRegionBased(this.Input, database, this.Overlap, this.Goroutine)
		if err != nil {
			return err
		}
		annoResults = append(annoResults, annoResult)
	}
	annoResult = anno.MergeAnnoResults(annoResults)
	vcfHeader.Infos = annoResult.VcfHeaderInfo
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	vcfWriter, err := vcfgo.NewWriter(writer, vcfHeader)
	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
		pk := (&pkg.Variant{Variant: *variant}).PK()
		for key, val := range annoResult.AnnoInfos[pk] {
			if val != "" && val != "." {
				err = variant.Info().Set(key, val)
				if err != nil {
					return err
				}
			}
		}
		vcfWriter.WriteVariant(variant)
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
			param.Goroutine, _ = cmd.Flags().GetInt("goroutine")
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
	cmd.Flags().IntP("goroutine", "c", 10000, "Parameter Goroutine Numbers")
	return cmd
}
