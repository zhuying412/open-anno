package pre

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"os"
	"path"
	"strings"

	"github.com/brentp/vcfgo"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type PreClinvarParam struct {
	Input  string `validate:"required,pathexists"`
	Output string `validate:"required"`
}

func (this PreClinvarParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this PreClinvarParam) Run() error {
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
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprint(writer, "#Chr\tStart\tEnd\tRef\tAlt\tCLNREVSTAT\tCLNSIG\tCLNSIGCONF\tCLNDN\n")
	for {
		row := vcfReader.Read()
		if row == nil {
			break
		}
		for _, alt := range row.Alt() {
			chrom, start, end, ref, alt := pkg.VCFtoAV(row.Chrom(), int(row.Pos), row.Ref(), alt)
			chrom = "chr" + chrom
			clnrevstat, err := row.Info().Get("CLNREVSTAT")
			if err != nil {
				clnrevstat = "."
			} else {
				if clnrevstats, ok := clnrevstat.([]string); ok {
					clnrevstat = strings.Join(clnrevstats, ",")
				}
			}
			clnsig, err := row.Info().Get("CLNSIG")
			if err != nil {
				clnsig = "."
			}
			clnsigconf, err := row.Info().Get("CLNSIGCONF")
			if err != nil {
				clnsigconf = "."
			} else {
				clnsigconf = strings.ReplaceAll(clnsigconf.(string), "|", ",")
			}
			clndn, err := row.Info().Get("CLNDN")
			if err != nil {
				clndn = "."
			} else {
				if clndns, ok := clndn.([]string); ok {
					clndn = strings.Join(clndns, ",")
				}
			}
			fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", chrom, start, end, ref, alt, clnrevstat, clnsig, clnsigconf, clndn)
		}
	}
	return err
}

func NewPreClinvarCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "cln",
		Short: "Prepare Clinvar Base on NCBI Clinvar VCF",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreClinvarParam
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
	cmd.Flags().StringP("input", "i", "", "Input Clinvar VCF File")
	cmd.Flags().StringP("output", "o", "", "Output File")
	return cmd
}
