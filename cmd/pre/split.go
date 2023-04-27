package pre

import (
	"fmt"
	"io"
	"log"
	"open-anno/pkg"
	"os"
	"strconv"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type PreSplitVCFParam struct {
	Input  string `validate:"required,pathexists"`
	Outdir string `validate:"required"`
}

func (this PreSplitVCFParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return os.MkdirAll(this.Outdir, 0666)
}

func (this PreSplitVCFParam) Run() error {
	reader, err := pkg.NewIOReader(this.Input)
	if err != nil {
		return err
	}
	defer reader.Close()
	writerMap := make(map[string]io.WriteCloser)
	scanner := pkg.NewIOScanner(reader)
	header := make([]string, 0)
	for scanner.Scan() {
		text := scanner.Text()
		if strings.HasPrefix(text, "#") {
			header = append(header, text)
		} else {
			row := strings.Split(text, "\t")
			pos, err := strconv.Atoi(row[1])
			if err != nil {
				return err
			}
			key := fmt.Sprintf("%s.%d", row[0], pos/pkg.FilterBasedBucketSize)
			if _, ok := writerMap[key]; !ok {
				writer, err := pkg.NewIOWriter(fmt.Sprintf("%s/%s.vcf", this.Outdir, key))
				if err != nil {
					return err
				}
				defer writer.Close()
				fmt.Fprintf(writer, "%s\n", strings.Join(header, "\n"))
				writerMap[key] = writer
			}
			fmt.Fprintf(writerMap[key], "%s\n", text)
		}
	}
	return err
}

func NewSplitVCFCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "splitvcf",
		Short: "Prepare Split VCF",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreSplitVCFParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.Outdir, _ = cmd.Flags().GetString("outdir")
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
	cmd.Flags().StringP("input", "i", "", "Input VCF File")
	cmd.Flags().StringP("outdir", "o", "", "Output Directory")
	return cmd
}
