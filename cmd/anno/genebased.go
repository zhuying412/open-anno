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

	"github.com/brentp/faidx"
	"github.com/spf13/cobra"
)

func initGeneBasedData(avinput string, dbPath string, dbName string, builder string) (map[string]variant.Variants, gene.Transcripts, gene.TransIndexes) {
	// builder
	gene.SetGenome(builder)
	// refgene
	refgeneFile := path.Join(dbPath, builder, dbName, "refgene.txt")
	log.Printf("Read Refgene: %s ...", refgeneFile)
	transcripts, err := gene.ReadRefgene(refgeneFile)
	if err != nil {
		log.Fatal(err)
	}
	// index
	indexFile := path.Join(dbPath, builder, dbName, "refgene.idx")
	log.Printf("Read Refgene Index: %s ...", indexFile)
	transIndexes, err := gene.ReadTransIndexs(indexFile)
	if err != nil {
		log.Fatal(err)
	}
	// snv
	log.Printf("Read avinput: %s ...", avinput)
	snvMap, err := variant.ReadAvinput(avinput)
	if err != nil {
		log.Fatal(err)
	}
	return snvMap, transcripts, transIndexes
}

func InitWriter(outfile string) *os.File {
	if outfile == "-" {
		return os.Stdout
	}
	writer, err := os.Create(outfile)
	if err != nil {
		log.Fatal(err)
	}
	return writer
}

func AnnoSnvGeneBased(avinput string, dbPath string, dbName string, builder string, outfile string, aashort bool) {
	snvMap, refgenes, refgeneIndexes := initGeneBasedData(avinput, dbPath, dbName, builder)
	// mrna
	mrnaFile := path.Join(dbPath, builder, dbName, "mRNA.fa")
	log.Printf("Read mRNA: %s ...", mrnaFile)
	fai, err := faidx.New(mrnaFile)
	if err != nil {
		log.Fatal(err)
	}
	// anno
	writer := InitWriter(outfile)
	fmt.Fprint(writer, "Chr\tStart\tEnd\tRef\tAlt\tGene\tGeneID\tEvent\tRegion\tDetail\n")
	for chrom, snvs := range snvMap {
		log.Printf("Start run annotate chr%s ...", chrom)
		transcripts := refgenes.FilterChrom(chrom, fai)
		transIndexes := refgeneIndexes.FilterChrom(chrom)
		if err != nil {
			log.Fatal(err)
		}
		genebased.AnnoSnvs(snvs, transcripts, transIndexes, aashort, writer)
	}
}

func AnnoCnvGeneBased(avinput string, dbPath string, dbName string, builder string, outfile string) {
	snvMap, refgenes, refgeneIndexes := initGeneBasedData(avinput, dbPath, dbName, builder)
	// anno
	writer := InitWriter(outfile)
	fmt.Fprint(writer, "Chr\tStart\tEnd\tRef\tAlt\tRegion\n")
	for chrom, cnvs := range snvMap {
		log.Printf("Filter GeneBased DB by %s ...", chrom)
		transcripts := refgenes.FilterChrom(chrom, &faidx.Faidx{})
		transIndexes := refgeneIndexes.FilterChrom(chrom)
		genebased.AnnoCnvs(cnvs, transcripts, transIndexes, writer)
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
