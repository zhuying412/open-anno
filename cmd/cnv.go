package cmd

import (
	"github.com/spf13/cobra"
	"grandanno/cnv"
	"log"
)

func NewCnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "cnv",
		Short: "Annotate CNV",
		Run: func(cmd *cobra.Command, args []string) {
			log.Println(cmd.Flag("database").Value.String())
			if cmd.Flag("database").Value.String() == "" ||
				cmd.Flag("input").Value.String() == "" ||
				cmd.Flag("output").Value.String() == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				cnv.AnnoCnv(
					cmd.Flag("database").Value.String(),
					cmd.Flag("genome").Value.String(),
					cmd.Flag("input").Value.String(),
					cmd.Flag("output").Value.String(),
				)
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "SNV input file")
	cmd.Flags().StringP("output", "o", "", "Annotated ouput file")
	cmd.Flags().StringP("database", "d", "", "Database Path")
	cmd.Flags().StringP("genome", "g", "hg19", "Genome Version")
	return cmd
}
