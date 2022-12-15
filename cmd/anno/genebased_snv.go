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
	Input    string `validate:"required,pathexists"`
	GenePred string `validate:"required,pathexists"`
	MRNA     string `validate:"required,pathexists"`
	RefDict  string `validate:"required,pathexists"`
	Gene     string `validate:"required,pathexists"`
	Output   string `validate:"required"`
	DBName   string `validate:"required"`
	AAshort  bool
	Exon     bool
}

func (this AnnoGBSnvParam) TransIndex() string {
	return this.GenePred + ".idx"
}

func (this AnnoGBSnvParam) Valid() error {
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
	log.Printf("Read TransIndex: %s ...", this.TransIndex())
	transIndexes, err := pkg.ReadTransIndexes(this.TransIndex())
	if err != nil {
		return err
	}
	// 读取mRNA输入文件
	log.Printf("Read mRNA: %s ...", this.MRNA)
	mrna, err := faidx.New(this.MRNA)
	if err != nil {
		return err
	}
	// 读取RefDict输入文件
	log.Printf("Read RefDict : %s ...", this.RefDict)
	genome, err := anno.ReadGenomeDict(this.RefDict)
	if err != nil {
		return err
	}
	// 读取gene输入文件
	log.Printf("Read Gene: %s ...", this.Gene)
	geneSymbolToID, err := anno.ReadGene(this.Gene)
	if err != nil {
		return err
	}
	return gene.AnnoSnvs(snvMap, this.Output, this.DBName, gpes, transIndexes, mrna, genome, geneSymbolToID, this.AAshort)
}

func NewAnnoGBSnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "snv",
		Short: "Annotate GeneBased of SNV",
		Run: func(cmd *cobra.Command, args []string) {
			var param AnnoGBSnvParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.MRNA, _ = cmd.Flags().GetString("mrna")
			param.RefDict, _ = cmd.Flags().GetString("refdict")
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
	cmd.Flags().StringP("genepred", "p", "", "Input GenePred File")
	cmd.Flags().StringP("mrna", "m", "", "Input mRNA File")
	cmd.Flags().StringP("refdict", "r", "", "Input Reference Dict File")
	cmd.Flags().StringP("gene", "g", "", "Input Gene Symbol To ID File")
	cmd.Flags().StringP("dbname", "d", "", "Parameter Database Name")
	cmd.Flags().StringP("output", "o", "", "AnnoOutput File")
	cmd.Flags().BoolP("aashort", "a", false, "Parameter Is AA Short")
	cmd.Flags().BoolP("exon", "e", false, "Parameter Is Exon")
	return cmd
}
