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

type AnnoFBParam struct {
	Input    string `validate:"required,pathexists"`
	Database string `validate:"required,pathexists"`
	Output   string `validate:"required"`
}

func (this AnnoFBParam) DBIdx() string {
	return this.Database + ".idx"
}

func (this AnnoFBParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this AnnoFBParam) Run() error {
	return db.AnnoFilterBased(this.Input, this.Database, this.DBIdx(), this.Output)
}

func NewAnnoFBCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "fb",
		Short: "Annotate FilterBased",
		Run: func(cmd *cobra.Command, args []string) {
			var param AnnoFBParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.Database, _ = cmd.Flags().GetString("database")
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
	cmd.Flags().StringP("database", "d", "", "Input Database File")
	cmd.Flags().StringP("output", "o", "", "AnnoOutput File")
	return cmd
}
