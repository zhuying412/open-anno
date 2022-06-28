package tools

import (
	"log"
	"open-anno/pkg/io"

	"github.com/spf13/cobra"
)

type Av2BedParam struct {
	Av2AnyParam
}

func (this Av2BedParam) Run() error {
	variants, err := this.Variants()
	if err != nil {
		return err
	}
	beds := make(io.BEDs, len(variants))
	for i, row := range variants {
		beds[i] = io.BED{Chrom: row.Chrom, Start: row.Start, End: row.End, Name: row.Alt}
	}
	return io.WriteBEDs(this.Ouput, beds)
}

func NewAV2BEDCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "av2bed",
		Short: "Convert CNV AVINPUT to BED",
		Run: func(cmd *cobra.Command, args []string) {
			var param Av2BedParam
			param.Input, _ = cmd.Flags().GetString("bed")
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
	cmd.Flags().StringP("avinput", "i", "", "AVINPUT File")
	cmd.Flags().StringP("bed", "o", "", "BED File")
	return cmd
}
