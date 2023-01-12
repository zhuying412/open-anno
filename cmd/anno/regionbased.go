package anno

import (
	"log"
	"open-anno/anno/db"
	"open-anno/pkg"
	"os"
	"path"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type AnnoRBParam struct {
	Input    string  `validate:"required,pathexists"`
	Database string  `validate:"required,pathexists"`
	Output   string  `validate:"required"`
	Overlap  float64 `validate:"required"`
}

func (this AnnoRBParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this AnnoRBParam) Run() error {
	return db.AnnoRegionBased(this.Input, this.Database, this.Output, this.Overlap)
}

func NewAnnoRBCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "rb",
		Short: "Annotate RegionBased",
		Run: func(cmd *cobra.Command, args []string) {
			var param AnnoRBParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.Database, _ = cmd.Flags().GetString("database")
			param.Output, _ = cmd.Flags().GetString("output")
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
	cmd.Flags().StringP("database", "d", "", "Input Database File")
	cmd.Flags().StringP("output", "o", "", "AnnoOutput File")
	cmd.Flags().Float64P("overlap", "l", 0.7, "Parameter Overlap")
	return cmd
}
