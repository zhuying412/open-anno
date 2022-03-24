package anno

import (
	"fmt"
	"log"
	"open-anno/anno/genebased"
	"open-anno/pkg/gene"
	"open-anno/pkg/variant"
	"os"
	"path"
	"strings"

	"github.com/spf13/cobra"
)

func AnnoSnvGeneBased(avinput string, dbPath string, dbName string, builder string, outfile string, aashort bool) {
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
	fmt.Fprint(writer, "Chr\tStart\tEnd\tRef\tAlt\tGene\tGeneID\tEvent\tRegion\tDetail\n")
	for chrom, snvs := range snv_dict {
		transFile := path.Join(dbPath, builder, dbName, fmt.Sprintf("chr%s.json", chrom))
		indexFile := path.Join(dbPath, builder, dbName, fmt.Sprintf("chr%s.idx.json", chrom))
		log.Printf("Read Transcript: %s ...", transFile)
		transcripts, err := gene.ReadTransDB(transFile)
		if err != nil {
			log.Fatal(err)
		}
		log.Printf("Read Transcript Index: %s ...", indexFile)
		transIndexes, err := gene.ReadTransIndexDB(indexFile)
		if err != nil {
			log.Fatal(err)
		}
		log.Printf("Start run annotate ...")
		genebased.AnnoSnvs(snvs, transcripts, transIndexes, aashort, writer)
	}
}

func AnnoCnvGeneBased(avinput string, dbPath string, dbName string, builder string, outfile string) {
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
	fmt.Fprint(writer, "Chr\tStart\tEnd\tRef\tAlt\tRegion\n")
	for chrom, snvs := range snv_dict {
		transFile := path.Join(dbPath, builder, dbName, fmt.Sprintf("chr%s.json", chrom))
		indexFile := path.Join(dbPath, builder, dbName, fmt.Sprintf("chr%s.idx.json", chrom))
		log.Printf("Read Transcript: %s ...", transFile)
		transcripts, err := gene.ReadTransDB(transFile)
		if err != nil {
			log.Fatal(err)
		}
		log.Printf("Read Transcript Index: %s ...", indexFile)
		transIndexes, err := gene.ReadTransIndexDB(indexFile)
		if err != nil {
			log.Fatal(err)
		}
		log.Printf("Start run annotate ...")
		genebased.AnnoCnvs(snvs, transcripts, transIndexes, writer)
	}
}

func NewAnnoGeneBasedCmd(varType string) *cobra.Command {
	varType = strings.ToLower(varType)
	cmd := &cobra.Command{
		Use:   fmt.Sprintf("%sGB", strings.Title(varType)),
		Short: fmt.Sprintf("Annotate %s Genebased", varType),
		Run: func(cmd *cobra.Command, args []string) {
			avinput, _ := cmd.Flags().GetString("avinput")
			dbpath, _ := cmd.Flags().GetString("dbpath")
			dbname, _ := cmd.Flags().GetString("dbname")
			builder, _ := cmd.Flags().GetString("builder")
			outfile, _ := cmd.Flags().GetString("outfile")
			if avinput == "" || dbpath == "" || dbname == "" || builder == "" || outfile == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				if varType == "snv" {
					aashort, _ := cmd.Flags().GetBool("aashort")
					AnnoSnvGeneBased(avinput, dbpath, dbname, builder, outfile, aashort)
				} else if varType == "CNV" {
					AnnoCnvGeneBased(avinput, dbpath, dbname, builder, outfile)
				}
			}
		},
	}
	cmd.Flags().StringP("avinput", "i", "", "Annotated Variants Input File")
	cmd.Flags().StringP("outfile", "o", "-", "Output File, - for stdout")
	cmd.Flags().StringP("dbpath", "d", "", "Database Directory")
	cmd.Flags().StringP("dbname", "n", "refgene", "Database Builder")
	cmd.Flags().StringP("builder", "b", "hg19", "Database Builder")
	if varType == "snv" {
		cmd.Flags().BoolP("aashort", "s", true, "Database Builder")
	}
	return cmd
}
