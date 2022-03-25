package tools

import (
	"log"
	"open-anno/pkg/variant"
	"open-anno/tools"

	"github.com/brentp/faidx"
	"github.com/spf13/cobra"
)

func RunVCf2AV(vcf string, avinput string) {
	variants, err := tools.ReadVCF(vcf)
	if err != nil {
		log.Fatal(err)
	}
	err = variant.WriteAvinput(variants, avinput)
	if err != nil {
		log.Fatal(err)
	}
}

func RunAV2VCF(avinput string, vcf string, genome string) {
	fai, err := faidx.New(genome)
	if err != nil {
		log.Fatal(err)
	}
	vcfs, err := tools.ReadAV(avinput, fai)
	if err != nil {
		log.Fatal(err)
	}
	err = tools.WriteVCF(vcfs, vcf)
	if err != nil {
		log.Fatal(err)
	}
}

func NewVCf2AVCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "vcf2av",
		Short: "Convert VCF to AVINPUT",
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
		Use:   "vcf2av",
		Short: "Convert VCF to AVINPUT",
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
