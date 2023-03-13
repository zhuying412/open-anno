package anno

import (
	"log"
	"open-anno/anno"
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
	Overlap            float64  `validate:"required"`
}

func (this *AnnoSnvParam) Valid() error {
	this.GenePredIndex = this.GenePred + ".idx"
	this.GenomeIndex = this.Genome + ".fai"
	for _, db := range this.FilterBaseds {
		this.FilterBasedIndexes = append(this.FilterBasedIndexes, db+".idx")
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

func (this *AnnoSnvParam) RunAnnoGeneBase(snvs anno.Variants, annoChan chan map[string]map[string]any, vcfHeaderInfoChan chan []vcfgo.Info, errChan chan error) {
	// 读取GenePred输入文件
	log.Printf("Read GenePred: %s ...", this.GenePred)
	gpes, err := pkg.ReadGenePred(this.GenePred)
	if err != nil {
		errChan <- err
		return
	}
	// 读取Transcript Index输入文件
	log.Printf("Read GenePred Index: %s ...", this.GenePredIndex)
	transIndexes, err := pkg.ReadTransIndexes(this.GenePredIndex)
	if err != nil {
		errChan <- err
		return
	}
	// 读取Genome输入文件
	log.Printf("Read Reference Fasta: %s ...", this.Genome)
	genome, err := faidx.New(this.Genome)
	if err != nil {
		errChan <- err
		return
	}
	// 读取gene输入文件
	log.Printf("Read Gene: %s ...", this.Gene)
	geneSymbolToID, err := anno.ReadGene(this.Gene)
	if err != nil {
		errChan <- err
		return
	}
	annoInfos, err := gene.AnnoSnvs(snvs, gpes, transIndexes, genome, geneSymbolToID, this.AAshort)
	vcfHeaderInfoChan <- []vcfgo.Info{
		{
			Id:          "GENE",
			Description: "Gene Symbol",
			Number:      ".",
			Type:        "String",
		},
		{
			Id:          "GENE_ID",
			Description: "Gene Entrez ID",
			Number:      ".",
			Type:        "String",
		},
		{
			Id:          "REGION",
			Description: "Region in gene, eg: exonic, intronic, UTR3, UTR5",
			Number:      ".",
			Type:        "String",
		},
		{
			Id:          "EVENT",
			Description: "Variant Event, eg: missense, nonsense, splicing",
			Number:      ".",
			Type:        "String",
		},
		{
			Id:          "DETAIL",
			Description: "Gene detail, FORMAT=Gene:Transcript:Exon:NA_CHANGE:AA_CHANGE",
			Number:      ".",
		},
	}
	annoChan <- annoInfos
	errChan <- err
}

func (this *AnnoSnvParam) RunAnnoFilterBase(snvs []anno.SNV, database string, annoChan chan map[string]map[string]any, vcfHeaderInfoChan chan []vcfgo.Info, errChan chan error) {
	annoInfos, headerInfos, err := db.AnnoFilterBased(snvs, database)
	annoChan <- annoInfos
	vcfHeaderInfoChan <- headerInfos
	errChan <- err
}

func (this AnnoSnvParam) Run() error {
	// 读取变异输入文件
	log.Printf("Read AnnoInput: %s ...", this.Input)
	snvs, vcfHeader, err := anno.ReadVCF(this.Input)
	if err != nil {
		return err
	}
	// 构造 error channel
	errChan := make(chan error, len(this.FilterBaseds)+1)
	annoInfosChan := make(chan map[string]map[string]any, len(this.FilterBaseds)+1)
	vcfHeaderInfoChan := make(chan []vcfgo.Info, len(this.FilterBaseds)+1)
	go this.RunAnnoGeneBase(snvs, annoInfosChan, vcfHeaderInfoChan, errChan)
	// FilterBased 注释
	for i, database := range this.FilterBaseds {
		go this.RunAnnoFilterBase(snvs, database, annoInfosChan, vcfHeaderInfoChan, errChan)
	}
	// 错误处理
	annoInfos := make(map[string]map[string]any)
	for i := 0; i < len(this.FilterBaseds); i++ {
		err := <-errChan
		if err != nil {
			return err
		}
		for pk, infos := range <-annoInfosChan {
			for key, val := range infos {
				_, ok := annoInfos[pk]
				if !ok {
					annoInfos[pk] = map[string]any{}
				}
				annoInfos[pk][key] = val
			}

		}
		for key, info := range <-vcfHeaderInfoChan {
			vcfHeader.Infos[key] = info
		}
	}
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()

	vcfWriter, err := vcfgo.NewWriter(writer, vcfHeader)
	for _, snv := range snvs {
		vcfVariant := snv.(anno.SNV).Variant
		for key, val := range annoInfos[snv.AnnoVariant().PK()] {
			err = vcfVariant.Info().Set(key, val)
			if err != nil {
				return err
			}
		}
		vcfWriter.WriteVariant(&vcfVariant)
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
	return cmd
}
