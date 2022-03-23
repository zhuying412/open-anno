package anno

import (
	"fmt"
	"log"
	"open-anno/anno/database"
	"open-anno/pkg/variant"
	"os"
	"path"

	"github.com/spf13/cobra"
)

func AnnoRegionBased(avinput string, dbPath string, dbName string, builder string, overlap float64, outfile string) {
	log.Printf("Read avinput: %s ...", avinput)
	snv_dict, err := variant.ReadAvinput(avinput)
	if err != nil {
		log.Fatal(err)
	}
	var writer *os.File
	if outfile == "-" {
		writer = os.Stdout
	} else {
		writer, err = os.Create(outfile)
		if err != nil {
			log.Fatal(err)
		}
	}
	for chrom, snvs := range snv_dict {
		dbfile := path.Join(dbPath, builder, dbName, fmt.Sprintf("chr%s.txt", chrom))
		database.AnnoRegionBased(snvs, dbfile, overlap, writer)
	}
}

func NewRegionBasedCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "RB",
		Short: "Annotate Regionbased",
		Run: func(cmd *cobra.Command, args []string) {
			avinput, _ := cmd.Flags().GetString("avinput")
			dbpath, _ := cmd.Flags().GetString("dbpath")
			dbname, _ := cmd.Flags().GetString("dbname")
			builder, _ := cmd.Flags().GetString("builder")
			outfile, _ := cmd.Flags().GetString("outfile")
			overlap, _ := cmd.Flags().GetFloat64("overlap")
			if avinput == "" || dbpath == "" || dbname == "" || builder == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				AnnoRegionBased(avinput, dbpath, dbname, builder, overlap, outfile)
			}
		},
	}
	cmd.Flags().StringP("avinput", "i", "", "Annotated Variants Input File")
	cmd.Flags().StringP("outfile", "o", "-", "Output File, - for stdout")
	cmd.Flags().StringP("dbpath", "d", "", "Database Directory")
	cmd.Flags().StringP("dbname", "n", "", "Database name")
	cmd.Flags().StringP("builder", "b", "hg19", "Database Builder")
	cmd.Flags().Float64P("overlap", "p", 0.75, "Overlap percent of snv")
	return cmd
}
