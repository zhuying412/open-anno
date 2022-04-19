package tools

import (
	"log"
	"open-anno/tools"

	"github.com/spf13/cobra"
)

func RunAV2BED(avinput string, bed string) {

	beds, err := tools.ReadCnvAV(avinput)
	if err != nil {
		log.Fatal(err)
	}
	err = tools.WriteBED(beds, bed)
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
