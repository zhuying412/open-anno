package tools

import (
	"log"
	"open-anno/pkg/io"

	"github.com/brentp/faidx"
	"github.com/spf13/cobra"
)

type Vcf2AvParam struct {
	Any2AnyParam
}

func (this Vcf2AvParam) Run() error {
	vcfs, err := io.ReadVCFs(this.Input)
	if err != nil {
		return err
	}
	variants := make(io.Variants, len(vcfs))
	for i, vcf := range vcfs {
		variants[i] = vcf.Variant()
	}
	return io.WriteVariants(this.Ouput, variants)
}

func NewVCf2AVCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "vcf2av",
		Short: "Convert SNV VCF to AVINPUT",
		Run: func(cmd *cobra.Command, args []string) {
			var param Vcf2AvParam
			param.Input, _ = cmd.Flags().GetString("vcf")
			param.Ouput, _ = cmd.Flags().GetString("avinput")
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
	cmd.Flags().StringP("vcf", "i", "", "VCF File")
	cmd.Flags().StringP("avinput", "o", "", "AVINPUT File")
	return cmd
}

type Av2VcfParam struct {
	Av2AnyParam
	Genome string `validate:"required"`
}

func (this Av2VcfParam) Run() error {
	fai, err := faidx.New(this.Genome)
	if err != nil {
		return err
	}
	variants, err := this.Variants()
	if err != nil {
		return err
	}
	vcfs := make(io.VCFs, len(variants))
	for i, row := range variants {
		vcfs[i], err = row.VCF(fai)
		return err
	}
	return io.WriteVCFs(this.Ouput, vcfs)
}

func NewAV2VCFCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "av2vcf",
		Short: "Convert SNV AVINPUT to VCF",
		Run: func(cmd *cobra.Command, args []string) {
			var param Av2VcfParam
			param.Input, _ = cmd.Flags().GetString("avinput")
			param.Ouput, _ = cmd.Flags().GetString("vcf")
			param.Genome, _ = cmd.Flags().GetString("genome")
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
	cmd.Flags().StringP("vcf", "i", "", "VCF File")
	cmd.Flags().StringP("avinput", "o", "", "AVINPUT File")
	cmd.Flags().StringP("genome", "g", "", "Genome Fasta File")
	return cmd
}
