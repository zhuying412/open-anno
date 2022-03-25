package pre

import (
	"fmt"
	"io"
	"log"
	"open-anno/pkg/gene"
	"os"
	"path"
	"strings"

	"github.com/spf13/cobra"
)

func PreGeneBased(refgene string, fasta string, builder string, indexStep int, outdir string) {
	log.Println("Init parameters ...")
	gene.SetGenome(builder)
	if _, err := os.Stat(outdir); os.IsNotExist(err) {
		err := os.MkdirAll(outdir, os.ModePerm)
		if err != nil {
			log.Fatal(err)
		}
	}
	outRefgene := path.Join(outdir, "refgene.txt")
	log.Printf("Copy refgene to %s ...", outRefgene)
	refgeneReader, err := os.Open(refgene)
	if err != nil {
		log.Fatal(err)
	}
	refgeneWriter, err := os.Create(outRefgene)
	if err != nil {
		log.Fatal(err)
	}
	defer refgeneReader.Close()
	defer refgeneWriter.Close()
	if _, err := io.Copy(refgeneWriter, refgeneReader); err != nil {
		log.Fatal(err)
	}
	log.Printf("Read refgene: %s ...", outRefgene)
	transcripts, err := gene.ReadRefgene(outRefgene)
	if err != nil {
		log.Fatal(err)
	}
	// log.Printf("Read genome: %s ...", fasta)
	// fai, err := faidx.New(fasta)
	// if err != nil {
	// 	log.Fatal(err)
	// }
	// outmRNA := path.Join(outdir, "mRNA.fa")
	// log.Printf("Write mRNA: %s ...", outmRNA)
	// transWriter, err := os.Create(outmRNA)
	// if err != err {
	// 	log.Fatal(err)
	// }
	// defer transWriter.Close()
	// for _, trans := range transcripts {
	// 	sequence, err := fai.Get(trans.Chrom, trans.TxStart-1, trans.TxEnd)
	// 	if err != nil {
	// 		sequence, err = fai.Get("chr"+trans.Chrom, trans.TxStart-1, trans.TxEnd)
	// 	}
	// 	sequence = strings.ToUpper(sequence)
	// 	fmt.Fprintf(transWriter, ">%s:%s:%s\n%s\n", trans.Chrom, trans.Gene, trans.Name, sequence)
	// }
	outIndex := path.Join(outdir, "refgene.idx")
	log.Printf("Write Transcript Index: %s ...", outIndex)
	idxWriter, err := os.Create(outIndex)
	if err != err {
		log.Fatal(err)
	}
	defer idxWriter.Close()
	transIndexes := gene.NewTransIndexes(indexStep)
	for _, index := range transIndexes {
		index.SetTranscripts(transcripts)
		if len(index.Transcripts) > 0 {
			fmt.Fprintf(idxWriter, "%s\t%d\t%d\t%s\n", index.Chrom, index.Start, index.End, strings.Join(index.Transcripts, ","))
		}
	}
}

func NewPreGeneBasedCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "GB",
		Short: "Prepare Genebased database",
		Run: func(cmd *cobra.Command, args []string) {
			genome, _ := cmd.Flags().GetString("genome")
			refgene, _ := cmd.Flags().GetString("refgene")
			dbpath, _ := cmd.Flags().GetString("dbpath")
			builder, _ := cmd.Flags().GetString("builder")
			name, _ := cmd.Flags().GetString("name")
			step, _ := cmd.Flags().GetInt("step")
			if genome == "" || refgene == "" || dbpath == "" || builder == "" || name == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				PreGeneBased(refgene, genome, builder, step, path.Join(dbpath, builder, name))
			}
		},
	}
	cmd.Flags().StringP("genome", "g", "", "Reference Fasta File")
	cmd.Flags().StringP("refgene", "r", "", "RefGene File")
	cmd.Flags().StringP("dbpath", "d", "", "Database Directory")
	cmd.Flags().StringP("name", "n", "", "Database Name")
	cmd.Flags().StringP("builder", "b", "hg19", "Database Path")
	cmd.Flags().IntP("step", "L", 300000, "Transcript Index Step Length")
	return cmd
}
