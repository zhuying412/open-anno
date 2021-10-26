package cmd

import (
	"github.com/spf13/cobra"
	"grandanno/snv"
	"log"
)

func NewSnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "snv",
		Short: "Annotate SNV",
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
				InitViper(cmd.Flag("database").Value.String())
				snv.AnnoSnv(cmd.Flag("input").Value.String(), cmd.Flag("output").Value.String())
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "SNV input file")
	cmd.Flags().StringP("output", "o", "", "Annotated ouput file")
	cmd.Flags().StringP("database", "d", "", "Database Path")
	return cmd
}
