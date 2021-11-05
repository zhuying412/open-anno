package cmd

import (
	"github.com/spf13/cobra"
	"grandanno/db"
	"log"
)

func IndexDatabase(databseFile string, databseIndexFile string) {
	// database file
	log.Printf("read %s\n", databseFile)
	indexes := db.ReadDatabaseFile(databseFile)
	// database index file
	log.Printf("write %s\n", databseIndexFile)
	db.CreateDatabseIndexFile(indexes, databseIndexFile)
}

func NewIndexDatabaseCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "index",
		Short: "Index database",
		Run: func(cmd *cobra.Command, args []string) {
			if cmd.Flag("input").Value.String() == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				databseFile := cmd.Flag("input").Value.String()
				databaseIndexFile := databseFile + ".idx"
				IndexDatabase(databseFile, databaseIndexFile)
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "Input database file")
	return cmd
}
