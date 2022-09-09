package tools

import (
	"log"
	"open-anno/cmd/anno"
	"open-anno/cmd/pre"
	"open-anno/pkg/io"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type MergeParam struct {
	AnnoInput    string   `validate:"required,pathexists"`
	AnnoGBOutput string   `validate:"required,pathexists"`
	AnnoOutputs  []string `validate:"required,pathsexists"`
	Output       string   `validate:"required"`
}

func (this MergeParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathsexists", anno.CheckPathsExists)
	validate.RegisterValidation("pathexists", pre.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return nil
}

func (this MergeParam) Run() error {
	return io.MergeAnnoResult(this.Output, this.AnnoInput, this.AnnoGBOutput, this.AnnoOutputs...)
}

func NewMergeCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "merge",
		Short: "Merge Annotation",
		Run: func(cmd *cobra.Command, args []string) {
			var param MergeParam
			param.AnnoInput, _ = cmd.Flags().GetString("input")
			param.AnnoGBOutput, _ = cmd.Flags().GetString("gbanno")
			param.AnnoOutputs, _ = cmd.Flags().GetStringArray("annos")
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
	cmd.Flags().StringP("input", "i", "", "Annotated Variants Input File")
	cmd.Flags().StringP("gbanno", "g", "", "Annotate Result File of GeneBased")
	cmd.Flags().StringArrayP("annos", "d", []string{}, "FilterBased or RegionBased Annotation Files")
	cmd.Flags().StringP("output", "o", "", "Output Merged Annotation File")
	return cmd
}
