package cmd

import (
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"grandanno/db"
	"grandanno/seq"
	"log"
	"os"
	"path"
)

func PrepareDatabase() {
	db.InitChromArray()
	db.InitDBParam()
	// reference
	refenceFastaFile := path.Join(viper.GetString("db.path"), viper.GetString("db.reference"))
	log.Printf("read %s\n", refenceFastaFile)
	reference := seq.ReadFastaFile(refenceFastaFile)
	// refgene
	refgeneFile := path.Join(viper.GetString("db.path"), viper.GetString("db.refgene"))
	log.Printf("read %s\n", refgeneFile)
	refgenes := db.ReadRefgeneFile(refgeneFile)
	// mRNA
	mrnaFastaDir := path.Join(viper.GetString("db.path"), viper.GetString("db.mrna_directory"))
	log.Printf("init and write mRNA to %s\n", mrnaFastaDir)
	if _, err := os.Stat(mrnaFastaDir); os.IsNotExist(err) {
		if err = os.MkdirAll(mrnaFastaDir, os.ModePerm); err != nil {
			log.Panic(err)
		}
	}
	mrnaFasta := db.NewMrnaFastaMap(refgenes, reference)
	for chrom, fasta := range mrnaFasta {
		mrnaFastaFile := path.Join(mrnaFastaDir, "chr"+chrom+".fa")
		seq.CreateFastaFile(fasta, mrnaFastaFile)
	}
	// Reference Index
	refIndexFile := path.Join(viper.GetString("db.path"), viper.GetString("db.refgene_index"))
	log.Printf("init and write %s\n", refIndexFile)
	db.CreateRefIndexFile(refgenes, refIndexFile)
}

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
				PrepareDatabase()
			}
		},
	}
	cmd.Flags().StringP("database", "d", "", "Database Path")
	return cmd
}
