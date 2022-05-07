package anno

import (
	"fmt"
	"log"
	"open-anno/anno/database"
	"open-anno/pkg/io"
	"open-anno/pkg/seq"
	"os"
	"path"
	"strings"

	"github.com/spf13/cobra"
)

func RunAnnoRegionBased(avinput string, dbPath string, dbName string, builder string, overlap float64, outfile string) {
	// builder
	seq.SetGenome(builder)
	// snvs
	log.Printf("Read avinput: %s ...", avinput)
	snvs, err := io.ReadVariantMap(avinput)
	if err != nil {
		log.Fatal(err)
	}
	// anno
	var writer *os.File
	if outfile == "-" {
		writer = os.Stdout
	} else {
		writer, err = os.Create(outfile)
		if err != nil {
			log.Fatal(err)
		}
	}
	headerWrited := true
	for chrom, subSnvs := range snvs {
		dbFile := path.Join(dbPath, builder, dbName, fmt.Sprintf("chr%s.txt", chrom))
		if _, err = os.Stat(dbFile); os.IsNotExist(err) {
			continue
		}
		database.AnnoRegionBased(subSnvs, dbFile, overlap, headerWrited, writer)
		headerWrited = false
	}
}

func NewRegionBasedCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "rb",
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
				RunAnnoRegionBased(avinput, dbpath, dbname, strings.ToLower(builder), overlap, outfile)
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
