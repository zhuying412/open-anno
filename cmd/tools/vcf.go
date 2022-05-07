package tools

import (
	"log"
	"open-anno/pkg/io"

	"github.com/brentp/faidx"
	"github.com/spf13/cobra"
)

func RunVCf2AV(vcf string, avinput string) {
	vcfs, err := io.ReadVCFs(vcf)
	if err != nil {
		log.Fatal(err)
	}
	variants := make(io.Variants, len(vcfs))
	for i, vcf := range vcfs {
		variants[i] = vcf.Variant()
	}
	err = io.WriteVariants(avinput, variants)
	if err != nil {
		log.Fatal(err)
	}
}

func RunAV2VCF(avinput string, vcf string, genome string) {
	fai, err := faidx.New(genome)
	if err != nil {
		log.Fatal(err)
	}
	variants, err := io.ReadVariants(avinput)
	if err != nil {
		log.Fatal(err)
	}
	vcfs := make(io.VCFs, len(variants))
	for i, row := range variants {
		vcfs[i] = row.VCF(fai)
	}
	err = io.WriteVCFs(vcf, vcfs)
	if err != nil {
		log.Fatal(err)
	}
}

func NewVCf2AVCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "vcf2av",
		Short: "Convert SNV VCF to AVINPUT",
		Run: func(cmd *cobra.Command, args []string) {
			vcf, _ := cmd.Flags().GetString("vcf")
			avinput, _ := cmd.Flags().GetString("avinput")
			if vcf == "" || avinput == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				RunVCf2AV(vcf, avinput)
			}
		},
	}
	cmd.Flags().StringP("vcf", "i", "", "VCF File")
	cmd.Flags().StringP("avinput", "o", "", "AVINPUT File")
	return cmd
}

func NewAV2VCFCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "av2vcf",
		Short: "Convert SNV AVINPUT to VCF",
		Run: func(cmd *cobra.Command, args []string) {
			vcf, _ := cmd.Flags().GetString("vcf")
			avinput, _ := cmd.Flags().GetString("avinput")
			genome, _ := cmd.Flags().GetString("genome")
			if vcf == "" || avinput == "" || genome == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				RunAV2VCF(avinput, vcf, genome)
			}
		},
	}
	cmd.Flags().StringP("vcf", "i", "", "VCF File")
	cmd.Flags().StringP("avinput", "o", "", "AVINPUT File")
	cmd.Flags().StringP("genome", "g", "", "Genome Fasta File")
	return cmd
}
