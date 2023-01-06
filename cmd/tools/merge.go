package tools

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"os"
	"path"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type MergeParam struct {
	AnnoInput      string   `validate:"required,pathexists"`
	GeneAnnoOutput string   `validate:"required,pathexists"`
	DBAnnoOutputs  []string `validate:"required,pathsexists"`
	Output         string   `validate:"required"`
}

func (this MergeParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathsexists", pkg.CheckPathsExists)
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this MergeParam) ReadAnnoInput() (map[string]string, error) {
	variants := make(map[string]string)
	reader, err := pkg.NewIOReader(this.AnnoInput)
	if err != nil {
		return variants, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		variant := strings.Join(row[0:5], "\t")
		variants[variant] = row[5]
	}
	return variants, nil
}

func (this MergeParam) ReadDBAnnos() ([]map[string]string, []string, error) {
	results := make([]map[string]string, len(this.DBAnnoOutputs))
	headers := make([]string, len(this.DBAnnoOutputs))
	for i, dbannoFile := range this.DBAnnoOutputs {
		reader, err := pkg.NewIOReader(dbannoFile)
		if err != nil {
			reader.Close()
			return results, headers, err
		}
		scanner := pkg.NewIOScanner(reader)
		scanner.Scan()
		headers[i] = strings.Join(strings.Split(scanner.Text(), "\t")[5:], "\t")
		result := make(map[string]string)
		for scanner.Scan() {
			row := strings.Split(scanner.Text(), "\t")
			variant := strings.Join(row[0:5], "\t")
			text := strings.Join(row[5:], "\t")
			result[variant] = text
		}
		results[i] = result
		reader.Close()
	}
	return results, headers, nil
}

func (this MergeParam) Run() error {
	variants, err := this.ReadAnnoInput()
	results, headers, err := this.ReadDBAnnos()
	if err != nil {
		return err
	}

	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	reader, err := pkg.NewIOReader(this.GeneAnnoOutput)
	if err != nil {
		return err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	scanner.Scan()
	fmt.Fprintf(writer,
		"Chr\tStart\tEnd\tRef\tAlt\t%s\tOtherinfo\n",
		strings.Join(append(strings.Split(scanner.Text(), "\t")[5:], headers...), "\t"),
	)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		variant := strings.Join(row[0:5], "\t")
		fmt.Fprint(writer, scanner.Text())
		for i, result := range results {
			if text, ok := result[variant]; ok {
				fmt.Fprintf(writer, "\t%s", text)
			} else {
				dots := pkg.NewArr(len(strings.Split(headers[i], "\t")), ".")
				fmt.Fprintf(writer, "\t%s", strings.Join(dots, "\t"))
			}
		}
		if info, ok := variants[variant]; ok {
			fmt.Fprintf(writer, "\t%s", info)
		} else {
			fmt.Fprint(writer, "\t.")
		}
		fmt.Fprint(writer, "\n")
	}
	return nil
}

func NewMergeCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "merge",
		Short: "Merge Annotation",
		Run: func(cmd *cobra.Command, args []string) {
			var param MergeParam
			param.AnnoInput, _ = cmd.Flags().GetString("input")
			param.GeneAnnoOutput, _ = cmd.Flags().GetString("geneanno")
			param.DBAnnoOutputs, _ = cmd.Flags().GetStringArray("dbannos")
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
	cmd.Flags().StringP("input", "i", "", "Input Annotated Variants Input File")
	cmd.Flags().StringP("geneanno", "g", "", "Input Annotate Result File of GeneBased")
	cmd.Flags().StringArrayP("dbannos", "d", []string{}, "Input FilterBased or RegionBased Annotation Files")
	cmd.Flags().StringP("output", "o", "", "Output Merged Annotation File")
	return cmd
}
