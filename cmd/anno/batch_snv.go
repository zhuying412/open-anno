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

type AnnoBatchSnvParam struct {
	Input        string `validate:"required,pathexists"`
	GenePred     string `validate:"required,pathexists"`
	Genome       string `validate:"required,pathexists"`
	Gene         string `validate:"required,pathexists"`
	Output       string `validate:"required"`
	GBName       string `validate:"required"`
	AAshort      bool
	Exon         bool
	FilterBaseds []string `validate:"pathsexists"`
	RegionBaseds []string `validate:"pathsexists"`
	Overlap      float64  `validate:"required"`
	Clean        bool
}

func (this *AnnoBatchSnvParam) Valid() error {
	pkg.IS_EXON_REGION = this.Exon
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	validate.RegisterValidation("pathsexists", pkg.CheckPathsExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return os.MkdirAll(this.Outdir(), 0666)
}

func (this AnnoBatchSnvParam) Sample() string {
	return strings.Split(path.Base(this.Input), ".")[0]
}

func (this AnnoBatchSnvParam) Outdir() string {
	return path.Dir(this.Output)
}

func (this *AnnoBatchSnvParam) RunAnnoGB(outfile string, errChan chan error) {
	gbParam := AnnoGBSnvParam{
		Input:    this.Input,
		GenePred: this.GenePred,
		Genome:   this.Genome,
		Gene:     this.Gene,
		DBName:   this.GBName,
		AAshort:  this.AAshort,
		Exon:     this.Exon,
		Output:   outfile + "_trans",
	}
	err := gbParam.Valid()
	if err != nil {
		errChan <- err
		return
	}
	err = gbParam.Run()
	if err != nil {
		errChan <- err
		return
	}
	aggsParam := tools.AggsParam{
		Input:  gbParam.Output,
		Output: outfile,
	}
	err = aggsParam.Valid()
	if err != nil {
		errChan <- err
		return
	}
	err = aggsParam.Run()
	if err != nil {
		errChan <- err
		return
	}
	if this.Clean {
		err = os.Remove(gbParam.Output)
		if err != nil {
			errChan <- err
			return
		}
	}
	errChan <- nil
}

func (this *AnnoBatchSnvParam) RunAnnoFB(outfile, database string, errChan chan error) {
	fbParam := AnnoFBParam{
		Input:    this.Input,
		Database: database,
		Output:   outfile,
	}
	err := fbParam.Valid()
	if err != nil {
		errChan <- err
		return
	}
	err = fbParam.Run()
	if err != nil {
		errChan <- err
		return
	}
	errChan <- nil
}

func (this *AnnoBatchSnvParam) RunAnnoRB(outfile, database string, errChan chan error) {
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

func (this AnnoBatchSnvParam) Run() error {
	// 构造 error channel
	errChan := make(chan error, len(this.FilterBaseds)+len(this.RegionBaseds)+1)
	// GeneBased 注释
	gbOutfile := path.Join(this.Outdir(), this.Sample()+".anno_output."+path.Base(this.GenePred))
	go this.RunAnnoGB(gbOutfile, errChan)
	// FilterBased 注释
	fbOutfiles := make([]string, len(this.FilterBaseds))
	for i, database := range this.FilterBaseds {
		outfile := path.Join(this.Outdir(), this.Sample()+".anno_output."+path.Base(database))
		fbOutfiles[i] = outfile
		go this.RunAnnoFB(outfile, database, errChan)
	}
	// RegionBased 注释
	rbOutfiles := make([]string, len(this.RegionBaseds))
	for i, database := range this.RegionBaseds {
		outfile := path.Join(this.Outdir(), this.Sample()+".anno_output."+path.Base(database))
		rbOutfiles[i] = outfile
		go this.RunAnnoRB(outfile, database, errChan)
	}
	// 错误处理
	for i := 0; i < len(this.FilterBaseds)+len(this.RegionBaseds)+1; i++ {
		err := <-errChan
		if err != nil {
			return err
		}
	}
	// 合并结果
	mergeParam := tools.MergeParam{
		AnnoInput:    this.Input,
		GeneBased:    gbOutfile,
		FilterBaseds: fbOutfiles,
		RegionBaseds: rbOutfiles,
		Output:       this.Output,
	}
	err := mergeParam.Valid()
	if err != nil {
		return err
	}
	err = mergeParam.Run()
	if err != nil {
		return err
	}
	if this.Clean {
		outfiles := []string{gbOutfile}
		outfiles = append(outfiles, fbOutfiles...)
		outfiles = append(outfiles, rbOutfiles...)
		for _, outfile := range outfiles {
			err = os.Remove(outfile)
			if err != nil {
				return err
			}
		}

	}
	return nil
}

func NewAnnoBatchSnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "snv",
		Short: "Annotate Batch of SNV",
		Run: func(cmd *cobra.Command, args []string) {
			var param AnnoBatchSnvParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.Genome, _ = cmd.Flags().GetString("genome")
			param.Gene, _ = cmd.Flags().GetString("gene")
			param.GBName, _ = cmd.Flags().GetString("gbname")
			param.Output, _ = cmd.Flags().GetString("output")
			param.AAshort, _ = cmd.Flags().GetBool("aashort")
			param.Exon, _ = cmd.Flags().GetBool("exon")
			param.FilterBaseds, _ = cmd.Flags().GetStringArray("filterbaseds")
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
	cmd.Flags().StringP("genome", "G", "", "Input Genome Fasta File")
	cmd.Flags().StringP("gene", "g", "", "Input Gene Symbol To ID File")
	cmd.Flags().StringP("gbname", "n", "", "Parameter Database Name")
	cmd.Flags().StringP("output", "o", "", "AnnoOutput File")
	cmd.Flags().BoolP("aashort", "a", false, "Parameter Is AA Short")
	cmd.Flags().BoolP("exon", "e", false, "Parameter Is Exon")
	cmd.Flags().StringArrayP("filterbaseds", "f", []string{}, "Input FilterBased Database File")
	cmd.Flags().StringArrayP("regionbaseds", "r", []string{}, "Input RegionBased Database File")
	cmd.Flags().Float64P("overlap", "l", 0.7, "Parameter Database Name")
	cmd.Flags().BoolP("clean", "c", false, "Clean Temporary File")
	return cmd
}
