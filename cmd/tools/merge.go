package tools

import (
	"bytes"
	"fmt"
	"io"
	"log"
	"open-anno/pkg"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type AnnoResult struct {
	ID   string `json:"id"`
	Text string `json:"text"`
}

type AnnoResultScanner struct {
	pkg.IOScanner
	FieldNames []string
}

func (this AnnoResultScanner) Header() string {
	return strings.Join(this.FieldNames, "\t")
}

func (this AnnoResultScanner) FillDot() string {
	var buffer bytes.Buffer
	for i := range this.FieldNames {
		if i == 0 {
			buffer.WriteString(".")
		} else {
			buffer.WriteString("\t.")
		}
	}
	return buffer.String()
}

func NewAnnoResultScanner(reader io.ReadCloser) AnnoResultScanner {
	scanner := pkg.NewIOScanner(reader)
	scanner.Scan()
	fieldNames := strings.Split(scanner.Text(), "\t")[5:]
	return AnnoResultScanner{IOScanner: scanner, FieldNames: fieldNames}
}

type AnnoResultScanners []AnnoResultScanner

func (this AnnoResultScanners) Header() string {
	fieldNames := make([]string, 0)
	for _, scanner := range this {
		fieldNames = append(fieldNames, scanner.FieldNames...)
	}
	return strings.Join(fieldNames, "\t")
}

func (this AnnoResultScanners) Results() []map[string]string {
	results := make([]map[string]string, len(this))
	for i, scanner := range this {
		result := make(map[string]string)
		for scanner.Scan() {
			row := scanner.Row()
			result[row.ID] = row.Text
		}
		results[i] = result
	}
	return results
}

func (this AnnoResultScanner) Row() AnnoResult {
	fields := strings.Split(this.Text(), "\t")
	return AnnoResult{
		ID:   strings.Join(fields[0:5], ":"),
		Text: strings.Join(fields[5:], "\t"),
	}
}

func MergeAnnoResult(outfile, annoInput, annoGBOutput string, annoOuputs ...string) error {
	writer, err := pkg.NewIOWriter(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	scanners := make(AnnoResultScanners, len(annoOuputs)+1)
	for i, annoOutput := range append(annoOuputs, annoInput) {
		reader, err := pkg.NewIOReader(annoOutput)
		if err != nil {
			return err
		}
		defer reader.Close()
		scanner := NewAnnoResultScanner(reader)
		scanners[i] = scanner
	}
	results := scanners.Results()
	reader, err := pkg.NewIOReader(annoGBOutput)
	if err != nil {
		return err
	}
	defer reader.Close()
	scanner := NewAnnoResultScanner(reader)
	fmt.Fprintf(writer, "Chr\tStart\tEnd\tRef\tAlt\t%s\t%s\n", scanner.Header(), scanners.Header())
	for scanner.Scan() {
		row := scanner.Row()
		fmt.Fprint(writer, scanner.Text())
		for i, result := range results {
			if text, ok := result[row.ID]; ok {
				fmt.Fprintf(writer, "\t%s", text)
			} else {
				fmt.Fprintf(writer, "\t%s", scanners[i].FillDot())
			}
		}
		fmt.Fprint(writer, "\n")
	}
	return nil
}

type MergeParam struct {
	AnnoInput    string   `validate:"required,pathexists"`
	AnnoGBOutput string   `validate:"required,pathexists"`
	AnnoOutputs  []string `validate:"required,pathsexists"`
	Output       string   `validate:"required"`
}

func (this MergeParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathsexists", pkg.CheckPathsExists)
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return nil
}

func (this MergeParam) Run() error {
	return MergeAnnoResult(this.Output, this.AnnoInput, this.AnnoGBOutput, this.AnnoOutputs...)
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
