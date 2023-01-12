package anno

import (
	"log"
	"open-anno/cmd/tools"
	"open-anno/pkg"
	"os"
	"path"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type AnnoBatchCnvParam struct {
	Input        string   `validate:"required,pathexists"`
	GenePred     string   `validate:"required,pathexists"`
	Gene         string   `validate:"required,pathexists"`
	Output       string   `validate:"required"`
	GBName       string   `validate:"required"`
	RegionBaseds []string `validate:"pathsexists"`
	Overlap      float64  `validate:"required"`
	Clean        bool
}

func (this *AnnoBatchCnvParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	validate.RegisterValidation("pathsexists", pkg.CheckPathsExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return os.MkdirAll(this.Outdir(), 0666)
}

func (this AnnoBatchCnvParam) Sample() string {
	return strings.Split(path.Base(this.Input), ".")[0]
}

func (this AnnoBatchCnvParam) Outdir() string {
	return path.Dir(this.Output)
}

func (this *AnnoBatchCnvParam) RunAnnoGB(outfile string, errChan chan error) {
	gbParam := AnnoGBCnvParam{
		Input:    this.Input,
		GenePred: this.GenePred,
		Gene:     this.Gene,
		DBName:   this.GBName,
		Output:   outfile,
	}
	err := gbParam.Valid()
	if err != nil {
		errChan <- err
		return
	}
	errChan <- nil
}

func (this *AnnoBatchCnvParam) RunAnnoRB(outfile, database string, errChan chan error) {
	rbParam := AnnoRBParam{
		Input:    this.Input,
		Database: database,
		Output:   outfile,
		Overlap:  this.Overlap,
	}
	err := rbParam.Valid()
	if err != nil {
		errChan <- err
		return
	}
	err = rbParam.Run()
	if err != nil {
		errChan <- err
		return
	}
	errChan <- nil
}

func (this AnnoBatchCnvParam) Run() error {
	// 构造 error channel
	errChan := make(chan error, len(this.RegionBaseds)+1)
	// GeneBased 注释
	gbOutfile := path.Join(this.Outdir(), this.Sample()+".anno_output."+path.Base(this.GenePred))
	go this.RunAnnoGB(gbOutfile, errChan)
	// RegionBased 注释
	rbOutfiles := make([]string, len(this.RegionBaseds))
	for i, database := range this.RegionBaseds {
		outfile := path.Join(this.Outdir(), this.Sample()+".anno_output."+path.Base(database))
		rbOutfiles[i] = outfile
		go this.RunAnnoRB(outfile, database, errChan)
	}
	// 错误处理
	for i := 0; i < len(this.RegionBaseds)+1; i++ {
		err := <-errChan
		if err != nil {
			return err
		}
	}
	// 合并结果
	mergeParam := tools.MergeParam{
		AnnoInput:    this.Input,
		GeneBased:    gbOutfile,
		RegionBaseds: rbOutfiles,
		Output:       this.Output,
	}
	err := mergeParam.Valid()
	if err != nil {
		return err
	}
	if this.Clean {
		outfiles := append(rbOutfiles, gbOutfile)
		for _, outfile := range outfiles {
			err = os.Remove(outfile)
			if err != nil {
				return err
			}
		}

	}
	return mergeParam.Run()
}

func NewAnnoBatchCnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "cnv",
		Short: "Annotate Batch of CNV",
		Run: func(cmd *cobra.Command, args []string) {
			var param AnnoBatchCnvParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.Gene, _ = cmd.Flags().GetString("gene")
			param.GBName, _ = cmd.Flags().GetString("gbname")
			param.Output, _ = cmd.Flags().GetString("output")
			param.RegionBaseds, _ = cmd.Flags().GetStringArray("regionbaseds")
			param.Overlap, _ = cmd.Flags().GetFloat64("overlap")
			param.Clean, _ = cmd.Flags().GetBool("clean")
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
	cmd.Flags().StringP("gbname", "n", "", "Parameter Database Name")
	cmd.Flags().StringP("output", "o", "", "AnnoOutput File")
	cmd.Flags().StringArrayP("regionbaseds", "r", []string{}, "Input RegionBased Database File")
	cmd.Flags().Float64P("overlap", "l", 0.7, "Parameter Database Name")
	cmd.Flags().BoolP("clean", "c", false, "Clean Temporary File")
	return cmd
}
