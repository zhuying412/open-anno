package anno

import (
	"log"
	"open-anno/anno/db"
	"open-anno/anno/gene"
	"open-anno/pkg"
	"os"
	"path"

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
	Overlap            float64  `validate:"required"`
	Goroutine          int      `validate:"required"`
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

func (this *AnnoSnvParam) RunAnnoGeneBase(annoInfoChan chan gene.AnnoInfos, vcfHeaderInfoChan chan map[string]*vcfgo.Info, errChan chan error) {
	// 读取Genome输入文件
	log.Printf("Read Reference Fasta: %s ...", this.Genome)
	genome, err := faidx.New(this.Genome)
	if err != nil {
		annoInfoChan <- gene.AnnoInfos{}
		vcfHeaderInfoChan <- map[string]*vcfgo.Info{}
		errChan <- err
		return
	}
	// 读取gene输入文件
	log.Printf("Read Gene: %s ...", this.Gene)
	err = pkg.InitGeneSymbolToID(this.Gene)
	if err != nil {
		annoInfoChan <- gene.AnnoInfos{}
		vcfHeaderInfoChan <- map[string]*vcfgo.Info{}
		errChan <- err
		return
	}
	log.Printf("Annotate %s", this.GenePred)
	annoInfos, headerInfos, err := gene.AnnoSnvs(this.Input, this.GenePred, genome, this.AAshort, this.Goroutine)
	vcfHeaderInfoChan <- headerInfos
	annoInfoChan <- annoInfos
	errChan <- err
}

func (this *AnnoSnvParam) RunAnnoFilterBase(database string, annoInfoChan chan gene.AnnoInfos, vcfHeaderInfoChan chan map[string]*vcfgo.Info, errChan chan error) {
	annoInfos, headerInfos, err := db.AnnoFilterBased(this.Input, database, this.Goroutine)
	annoInfoChan <- annoInfos
	vcfHeaderInfoChan <- headerInfos
	errChan <- err
}

func (this *AnnoSnvParam) RunAnnoRegionBase(database string, annoInfoChan chan gene.AnnoInfos, vcfHeaderInfoChan chan map[string]*vcfgo.Info, errChan chan error) {
	annoInfos, headerInfos, err := db.AnnoRegionBased(this.Input, database, this.Overlap, this.Goroutine)
	annoInfoChan <- annoInfos
	vcfHeaderInfoChan <- headerInfos
	errChan <- err
}

func (this AnnoSnvParam) Run() error {
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
	annoSize := len(this.FilterBaseds) + len(this.RegionBaseds) + 1
	errChan := make(chan error, annoSize)
	annoInfosChan := make(chan gene.AnnoInfos, annoSize)
	vcfHeaderInfoChan := make(chan map[string]*vcfgo.Info, annoSize)
	go this.RunAnnoGeneBase(annoInfosChan, vcfHeaderInfoChan, errChan)
	// FilterBased 注释
	for _, database := range this.FilterBaseds {
		go this.RunAnnoFilterBase(database, annoInfosChan, vcfHeaderInfoChan, errChan)
	}
	// RegionBased 注释
	for _, database := range this.RegionBaseds {
		go this.RunAnnoRegionBase(database, annoInfosChan, vcfHeaderInfoChan, errChan)
	}
	// 错误处理
	annoInfos := make(map[string]map[string]any)
	for i := 0; i < annoSize; i++ {
		err := <-errChan
		if err != nil {
			return err
		}
		for pk, infos := range <-annoInfosChan {
			for key, val := range infos {
				if val != "." && val != "" {
					_, ok := annoInfos[pk]
					if !ok {
						annoInfos[pk] = map[string]any{}
					}
					annoInfos[pk][key] = val
				}
			}

		}
		for key, info := range <-vcfHeaderInfoChan {
			vcfHeader.Infos[key] = info
		}
	}
	close(errChan)
	close(annoInfosChan)
	close(vcfHeaderInfoChan)
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	vcfWriter, err := vcfgo.NewWriter(writer, vcfHeader)
	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
		pk := (&pkg.Variant{Variant: *variant}).PK()
		for key, val := range annoInfos[pk] {
			err = variant.Info().Set(key, val)
			if err != nil {
				return err
			}
		}
		vcfWriter.WriteVariant(variant)
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
	cmd.Flags().StringP("genome", "G", "", "Input Genome Fasta File")
	cmd.Flags().StringP("gene", "g", "", "Input Gene Symbol To ID File")
	cmd.Flags().StringP("output", "o", "", "AnnoOutput File")
	cmd.Flags().BoolP("aashort", "a", false, "Parameter Is AA Short")
	cmd.Flags().BoolP("exon", "e", false, "Parameter Is Exon")
	cmd.Flags().StringArrayP("filterbaseds", "f", []string{}, "Input FilterBased Database File")
	cmd.Flags().Float64P("overlap", "l", 0.7, "Parameter Database Name")
	cmd.Flags().IntP("goroutine", "c", 10000, "Parameter Goroutine Numbers")
	return cmd
}
