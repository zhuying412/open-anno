package tools

import (
	"log"
	"open-anno/cmd/pre"
	"open-anno/pkg/io"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type AggsParam struct {
	AnnoGBOutput string `validate:"required,pathexists"`
	DBname       string `validate:"required"`
	Output       string `validate:"required"`
}

func (this AggsParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pre.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return nil
}

func (this AggsParam) Run() error {
	return io.AggsSnvGBAnno(this.AnnoGBOutput, this.DBname, this.Output)
}

func NewAggsCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "aggs",
		Short: "Aggregate Snv GeneBased Annotation",
		Run: func(cmd *cobra.Command, args []string) {
			var param AggsParam
			param.AnnoGBOutput, _ = cmd.Flags().GetString("gbanno")
			param.DBname, _ = cmd.Flags().GetString("dbname")
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
	cmd.Flags().StringP("gbanno", "g", "", "Annotate Result File of GeneBased")
	cmd.Flags().StringP("dbname", "n", "", "DB Name of GeneBased")
	cmd.Flags().StringP("output", "o", "", "Output Aggregated Annotation File")
	return cmd
}
