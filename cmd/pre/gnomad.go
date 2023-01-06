package pre

import (
	"fmt"
	"io"
	"log"
	"open-anno/pkg"
	"os"
	"path"
	"strings"

	"github.com/brentp/vcfgo"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type PreGnomadParam struct {
	Input  string `validate:"required,pathexists"`
	Output string `validate:"required"`
}

func (this PreGnomadParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this PreGnomadParam) KeyGenerator(dbname string) []string {
	populations := []string{"", "oth", "ami", "sas", "fin", "eas", "amr", "afr", "asj", "nfe"}
	prefixes := []string{"AC", "AN", "AF", "nhomalt"}
	keys := make([]string, 0)
	for _, population := range populations {
		for _, prefix := range prefixes {
			key := prefix
			if population != "" {
				key += "_" + population
			}
			if dbname != "" {
				key = dbname + "_" + key
			}
			keys = append(keys, key)
		}
	}
	return keys

}

func (this PreGnomadParam) ProcessVCF(vcf string, keys []string, writer io.WriteCloser) error {
	reader, err := pkg.NewIOReader(this.Input)
	if err != nil {
		return err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return err
	}
	defer vcfReader.Close()
	for {
		row := vcfReader.Read()
		if row == nil {
			break
		}
		for _, alt := range row.Alt() {
			chrom, start, end, ref, alt := pkg.VCFtoAV(row.Chrom(), int(row.Pos), row.Ref(), alt)
			fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t", chrom, start, end, ref, alt)
			values := make([]interface{}, 0)
			for i, key := range keys {
				val, err := row.Info().Get(key)
				if err != nil {
					return err
				}
				values[i] = val
			}
			fmt.Fprintf(writer, "%s\n", "")
		}
	}
	return nil
}

func (this PreGnomadParam) Run() error {
	infoKeys := this.KeyGenerator("")
	mapKeys := this.KeyGenerator("gnomAD")
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprintf(writer, "#Chr\tStart\tEnd\tRef\tAlt\t%s\n", strings.Join(mapKeys, "\t"))
	vcfs := make([]string, 0)
	for _, vcf := range vcfs {
		err := this.ProcessVCF(vcf, infoKeys, writer)
		if err != nil {
			return err
		}
	}
	return nil
}

func NewPreGnomadCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "gnomad",
		Short: "Prepare GnomAD Frequency Base on GnomAD VCF",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreGnomadParam
			param.Input, _ = cmd.Flags().GetString("input")
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
	cmd.Flags().StringP("input", "i", "", "Input 1000Genomes VCF File")
	cmd.Flags().StringP("output", "o", "", "Output File")
	return cmd
}
