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
	AVinput string   `validate:"required,pathexists"`
	GBAnno  string   `validate:"required,pathexists"`
	DBAnnos []string `validate:"required,pathsexists"`
	Output  string   `validate:"required"`
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
	return io.MergeAnno(this.Output, this.AVinput, this.GBAnno, this.DBAnnos...)
}

func NewMergeCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "merge",
		Short: "Merge Annotation",
		Run: func(cmd *cobra.Command, args []string) {
			var param MergeParam
			param.AVinput, _ = cmd.Flags().GetString("avinput")
			param.GBAnno, _ = cmd.Flags().GetString("gbanno")
			param.DBAnnos, _ = cmd.Flags().GetStringArray("dbannos")
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
	cmd.Flags().StringP("avinput", "i", "", "Annotated Variants Input File")
	cmd.Flags().StringP("gbanno", "g", "", "GeneBased Annotation File")
	cmd.Flags().StringArrayP("dbannos", "d", []string{}, "FilterBased or RegionBased Annotation Files")
	cmd.Flags().StringP("output", "o", "", "Output Merged Annotation File")
	return cmd
}
