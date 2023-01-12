package anno

import (
	"log"
	"open-anno/anno"
	"open-anno/anno/gene"
	"open-anno/pkg"
	"os"
	"path"

	"github.com/brentp/faidx"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type AnnoGBSnvParam struct {
	Input         string `validate:"required,pathexists"`
	GenePred      string `validate:"required,pathexists"`
	GenePredIndex string `validate:"required,pathexists"`
	Genome        string `validate:"required,pathexists"`
	GenomeIndex   string `validate:"required,pathexists"`
	Gene          string `validate:"required,pathexists"`
	Output        string `validate:"required"`
	DBName        string `validate:"required"`
	AAshort       bool
	Exon          bool
}

func (this *AnnoGBSnvParam) Valid() error {
	this.GenePredIndex = this.GenePred + ".idx"
	this.GenomeIndex = this.Genome + ".fai"
	pkg.IS_EXON_REGION = this.Exon
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this AnnoGBSnvParam) Run() error {
	// 读取变异输入文件
	log.Printf("Read AnnoInput: %s ...", this.Input)
	snvMap, err := anno.ReadAnnoInput(this.Input)
	if err != nil {
		return err
	}
	// 读取GenePred输入文件
	log.Printf("Read GenePred: %s ...", this.GenePred)
	gpes, err := pkg.ReadGenePred(this.GenePred)
	if err != nil {
		return err
	}
	// 读取Transcript Index输入文件
	log.Printf("Read GenePred Index: %s ...", this.GenePredIndex)
	transIndexes, err := pkg.ReadTransIndexes(this.GenePredIndex)
	if err != nil {
		return err
	}
	// 读取Genome输入文件
	log.Printf("Read Reference Fasta: %s ...", this.Genome)
	genome, err := faidx.New(this.Genome)
	if err != nil {
		return err
	}
	// 读取gene输入文件
	log.Printf("Read Gene: %s ...", this.Gene)
	geneSymbolToID, err := anno.ReadGene(this.Gene)
	if err != nil {
		return err
	}
	return gene.AnnoSnvs(snvMap, this.Output, this.DBName, gpes, transIndexes, genome, geneSymbolToID, this.AAshort)
}

func NewAnnoGBSnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "snv",
		Short: "Annotate GeneBased of SNV",
		Run: func(cmd *cobra.Command, args []string) {
			var param AnnoGBSnvParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.Genome, _ = cmd.Flags().GetString("genome")
			param.Gene, _ = cmd.Flags().GetString("gene")
			param.DBName, _ = cmd.Flags().GetString("dbname")
			param.Output, _ = cmd.Flags().GetString("output")
			param.AAshort, _ = cmd.Flags().GetBool("aashort")
			param.Exon, _ = cmd.Flags().GetBool("exon")
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
	cmd.Flags().StringP("dbname", "n", "", "Parameter Database Name")
	cmd.Flags().StringP("output", "o", "", "AnnoOutput File")
	cmd.Flags().BoolP("aashort", "a", false, "Parameter Is AA Short")
	cmd.Flags().BoolP("exon", "e", false, "Parameter Is Exon")
	return cmd
}
