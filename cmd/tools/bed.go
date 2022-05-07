package tools

import (
	"log"
	"open-anno/pkg/io"

	"github.com/spf13/cobra"
)

func RunAV2BED(avinput string, bed string) {
	variants, err := io.ReadVariants(avinput)
	if err != nil {
		log.Fatal(err)
	}
	beds := make(io.BEDs, len(variants))
	for i, row := range variants {
		beds[i] = io.BED{Chrom: row.Chrom, Start: row.Start, End: row.End, Name: row.Alt}
	}
	err = io.WriteBEDs(bed, beds)
	if err != nil {
		log.Fatal(err)
	}
}

func NewAV2BEDCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "av2vcf",
		Short: "Convert CNV AVINPUT to BED",
		Run: func(cmd *cobra.Command, args []string) {
			bed, _ := cmd.Flags().GetString("bed")
			avinput, _ := cmd.Flags().GetString("avinput")
			if bed == "" || avinput == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				RunAV2BED(avinput, bed)
			}
		},
	}
	cmd.Flags().StringP("vcf", "i", "", "VCF File")
	cmd.Flags().StringP("avinput", "o", "", "AVINPUT File")
	return cmd
}
