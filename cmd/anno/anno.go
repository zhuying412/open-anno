package anno

import (
	"errors"
	"fmt"
	"log"
	"open-anno/cmd/pre"
	"open-anno/pkg/seq"
	"path"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

var (
	VType_SNV string = "SNV"
	VType_CNV string = "CNV"
)

var (
	DType_G string = "GeneBased"
	DType_F string = "FilterBased"
	DType_R string = "RegionBased"
)

type IAnnoParam interface {
	Bind(*pflag.FlagSet)
	Valid() error
	Run() error
}

type AnnoParam struct {
	Input     string `validate:"required,pathexists"`
	DBpath    string `validate:"required,pathexists"`
	DBname    string `validate:"required"`
	Builder   string `validate:"required"`
	OutPrefix string `validate:"required"`
}

func (this AnnoParam) Output() string {
	return fmt.Sprintf("%s.%s.anno.txt", this.OutPrefix, this.DBname)
}

func (this AnnoParam) DBFile() string {
	return path.Join(this.DBpath, this.Builder, this.DBname+".txt")
}

func (this AnnoParam) DBIndex() string {
	return path.Join(this.DBpath, this.Builder, this.DBname+".txt.idx")
}

func (this AnnoParam) GeneInfo() string {
	return path.Join(this.DBpath, this.Builder, this.DBname+".geneinfo.txt")
}

func (this *AnnoParam) Bind(flagSet *pflag.FlagSet) {
	this.Input, _ = flagSet.GetString("avinput")
	this.OutPrefix, _ = flagSet.GetString("outprefix")
	this.DBpath, _ = flagSet.GetString("dbpath")
	this.Builder, _ = flagSet.GetString("builder")
	this.DBname, _ = flagSet.GetString("dbname")
}

func (this AnnoParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pre.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	seq.SetGenome(this.Builder)
	return nil
}

func (this AnnoParam) Run() error {
	return errors.New("Not Impl")
}

func NewAnnoCmd(use, varType, dbType string) *cobra.Command {
	cmd := &cobra.Command{
		Use:   use,
		Short: fmt.Sprintf("Annotate %s for %s", dbType, varType),
		Run: func(cmd *cobra.Command, args []string) {
			var param IAnnoParam
			switch varType {
			case VType_SNV:
				switch dbType {
				case DType_G:
					param = new(AnnoSnvGBParam)
				case DType_F:
					param = new(AnnoSnvFBParam)
				case DType_R:
					param = new(AnnoSnvRBParam)
				}
			case VType_CNV:
				switch dbType {
				case DType_G:
					param = new(AnnoCnvGBParam)
				case DType_R:
					param = new(AnnoCnvRBParam)
				}
			}
			param.Bind(cmd.Flags())
			fmt.Println(param)
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
	cmd.Flags().StringP("outprefix", "o", "", "Output Prefix")
	cmd.Flags().StringP("dbpath", "d", "", "Database Directory")
	cmd.Flags().StringP("dbname", "n", "", "Database Names")
	cmd.Flags().StringP("builder", "b", "hg38", "Database Builder")
	if varType == VType_SNV && dbType == DType_G {
		cmd.Flags().BoolP("aashort", "s", false, "Database Builder")
		cmd.Flags().BoolP("exon", "e", false, "Output ExonOrder Instead of TypeOrder")
	}
	if dbType == DType_R {
		cmd.Flags().Float64P("overlap", "p", 0.75, "CNV Overlap Threshold")
	}
	return cmd
}
