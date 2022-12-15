package pre

import (
	"log"
	"open-anno/pkg"
	"os"
	"path"

	"github.com/brentp/faidx"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type PreMRNAParam struct {
	Reference string `validate:"required,pathexists"`
	GenePred  string `validate:"required,pathexists"`
	MRNA      string `validate:"required"`
}

func (this PreMRNAParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.MRNA)
	return os.MkdirAll(outdir, 0666)
}

func (this PreMRNAParam) Run() error {
	log.Printf("Read GenePred: %s ...", this.GenePred)
	gpes, err := pkg.ReadGenePred(this.GenePred)
	if err != err {
		return err
	}
	log.Printf("Read genome: %s ...", this.Reference)
	fai, err := faidx.New(this.Reference)
	if err != err {
		return err
	}
	return pkg.CreateMRNA(gpes, fai, this.MRNA)
}

func NewPreMRNACmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "mrna",
		Short: "Prepare mRNA FASTA File",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreMRNAParam
			param.Reference, _ = cmd.Flags().GetString("reference")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.MRNA, _ = cmd.Flags().GetString("mrna")
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
	cmd.Flags().StringP("reference", "r", "", "Input Reference Fasta File")
	cmd.Flags().StringP("genepred", "g", "", "Input RefGene File")
	cmd.Flags().StringP("mrna", "o", "", "Output mRNA File")
	return cmd
}
