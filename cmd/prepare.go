package cmd

import (
	"github.com/spf13/cobra"
	"grandanno/db"
	"log"
)

func NewPrepareCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "prepare",
		Short: "Prepare required database files",
		Run: func(cmd *cobra.Command, args []string) {
			log.Println(cmd.Flag("database").Value.String())
			if cmd.Flag("database").Value.String() == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				InitViper(cmd.Flag("database").Value.String())
				db.PrepareDatabase()
			}
		},
	}
	cmd.Flags().StringP("database", "d", "", "Database Path")
	return cmd
}
