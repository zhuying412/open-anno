package anno

import (
	"log"
	"open-anno/anno"
	"open-anno/anno/gene"
	"open-anno/pkg"
	"os"
	"path"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type AnnoGBCnvParam struct {
	Input         string `validate:"required,pathexists"`
	GenePred      string `validate:"required,pathexists"`
	GenePredIndex string `validate:"required,pathexists"`
	Gene          string `validate:"required,pathexists"`
	Output        string `validate:"required"`
	DBName        string `validate:"required"`
}

func (this *AnnoGBCnvParam) Valid() error {
	this.GenePredIndex = this.GenePred + ".idx"
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this AnnoGBCnvParam) Run() error {
	// 读取变异输入文件
	log.Printf("Read AnnoInput: %s ...", this.Input)
	cnvMap, err := anno.ReadAnnoInput(this.Input)
	if err != nil {
		return err
	}
	// 读取GenePred输入文件
	log.Printf("Read GenePred: %s ...", this.GenePred)
	gpes, err := pkg.ReadGenePred(this.GenePred)
	if err != nil {
		return err
	}
	// 读取GenePred Index输入文件
	log.Printf("Read GenePred Index: %s ...", this.GenePredIndex)
	transIndexes, err := pkg.ReadTransIndexes(this.GenePredIndex)
	if err != nil {
		return err
	}
	// 读取gene输入文件
	log.Printf("Read Gene: %s ...", this.Gene)
	geneSymbolToID, err := anno.ReadGene(this.Gene)
	if err != nil {
		return err
	}
	return gene.AnnoCnvs(cnvMap, gpes, transIndexes, geneSymbolToID, this.Output, this.DBName)
}

func NewAnnoGBCnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "cnv",
		Short: "Annotate GeneBased of CNV",
		Run: func(cmd *cobra.Command, args []string) {
			var param AnnoGBCnvParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.Gene, _ = cmd.Flags().GetString("gene")
			param.DBName, _ = cmd.Flags().GetString("dbname")
			param.Output, _ = cmd.Flags().GetString("output")
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
	cmd.Flags().StringP("dbname", "n", "", "Parameter Database Name")
	cmd.Flags().StringP("output", "o", "", "AnnoOutput File")
	return cmd
}
