package prepare

import (
	"OpenAnno/db/database"
	"github.com/spf13/cobra"
	"log"
)

func NewPrepareDatabaseCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "database",
		Short: "Prepare required database files",
		Run: func(cmd *cobra.Command, args []string) {
			input, _ := cmd.Flags().GetString("input")
			outdir, _ := cmd.Flags().GetString("outdir")
			if input == "" || outdir == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				database.Generate(input, outdir)
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "Input Database File")
	cmd.Flags().StringP("outdir", "o", "", "Output Database Directory")
	return cmd
}
