package anno

import (
	"log"
	"open-anno/cmd/pre"
	"open-anno/pkg/io"
	"os"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

func CheckPathsExists(fl validator.FieldLevel) bool {
	for _, path := range fl.Field().Interface().([]string) {
		_, err := os.Stat(path)
		if os.IsNotExist(err) {
			return false
		}
	}
	return true
}

type MergeParam struct {
	AnnoInput   string   `validate:"required,pathexists"`
	AnnoOutputs []string `validate:"required,pathsexists"`
	Output      string   `validate:"required"`
}

func (this MergeParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathsexists", CheckPathsExists)
	validate.RegisterValidation("pathexists", pre.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return nil
}

func (this MergeParam) Run() error {
	return io.MergeAnnoResult(this.Output, this.AnnoInput, this.AnnoOutputs...)
}

func NewMergeCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "merge",
		Short: "Merge Annotation",
		Run: func(cmd *cobra.Command, args []string) {
			var param MergeParam
			param.AnnoInput, _ = cmd.Flags().GetString("annoInput")
			param.AnnoOutputs, _ = cmd.Flags().GetStringArray("annoOutputs")
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
	cmd.Flags().StringP("annoInput", "i", "", "Annotated Variants Input File")
	cmd.Flags().StringArrayP("annoOutputs", "d", []string{}, "FilterBased or RegionBased Annotation Files")
	cmd.Flags().StringP("output", "o", "", "Output Merged Annotation File")
	return cmd
}
