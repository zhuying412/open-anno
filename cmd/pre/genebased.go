package pre

import (
	"fmt"
	"io"
	"log"
	"open-anno/pkg/gene"
	"open-anno/pkg/seq"
	"os"
	"os/exec"
	"path"
	"strings"

	"github.com/brentp/faidx"
	"github.com/spf13/cobra"
)

func RunPreGeneBased(refgene string, maneSelect string, ncbiGeneInfo string, fasta string, builder string, indexStep int, outdir string) {
	// param
	log.Println("Init parameters ...")
	gene.SetGenome(builder)
	if _, err := os.Stat(outdir); os.IsNotExist(err) {
		err := os.MkdirAll(outdir, os.ModePerm)
		if err != nil {
			log.Fatal(err)
		}
	}
	// entrez_id
	outGeneId := path.Join(outdir, "geneid.txt")
	log.Printf("Init gene symbol to entrez_id: %s", outGeneId)
	geneSymbolToId, err := gene.NewGeneSymbolToId(maneSelect, ncbiGeneInfo, refgene)
	if err != nil {
		log.Fatal(err)
	}
	writer, err := os.Create(outGeneId)
	if err != err {
		log.Fatal(err)
	}
	for symbol, entrezId := range geneSymbolToId {
		fmt.Fprintf(writer, "%s\t%s\n", symbol, entrezId)
	}
	os.Exit(0)
	// refgene
	outRefgene := path.Join(outdir, "refgene.txt")
	log.Printf("Copy refgene to %s ...", outRefgene)
	reader, err := os.Open(refgene)
	if err != nil {
		log.Fatal(err)
	}
	writer, err = os.Create(outRefgene)
	if err != nil {
		log.Fatal(err)
	}
	if _, err := io.Copy(writer, reader); err != nil {
		log.Fatal(err)
	}
	log.Printf("Read refgene: %s ...", outRefgene)
	transcripts, err := gene.ReadRefgene(outRefgene)
	if err != nil {
		log.Fatal(err)
	}
	// genome
	log.Printf("Read genome: %s ...", fasta)
	fai, err := faidx.New(fasta)
	if err != nil {
		log.Fatal(err)
	}
	// mRNA
	outmRNA := path.Join(outdir, "mRNA.fa")
	log.Printf("Write mRNA: %s ...", outmRNA)
	writer, err = os.Create(outmRNA)
	if err != err {
		log.Fatal(err)
	}
	for _, trans := range transcripts {
		sequence, err := seq.Fetch(fai, trans.Chrom, trans.TxStart-1, trans.TxEnd)
		if err != nil {
			log.Fatal(err)
		}
		sequence = strings.ToUpper(sequence)
		fmt.Fprintf(writer, ">%s:%s:%s\n%s\n", trans.Chrom, trans.Gene, trans.Name, sequence)
	}
	command := exec.Command("samtools", "faidx", outmRNA)
	err = command.Run()
	if err != nil {
		log.Fatal(err)
	}
	// index
	outIndex := path.Join(outdir, "refgene.idx")
	log.Printf("Write Transcript Index: %s ...", outIndex)
	writer, err = os.Create(outIndex)
	if err != err {
		log.Fatal(err)
	}
	transIndexes := gene.NewTransIndexes(indexStep)
	for _, index := range transIndexes {
		index.SetTranscripts(transcripts)
		if len(index.Transcripts) > 0 {
			fmt.Fprintf(writer, "%s\t%d\t%d\t%s\n", index.Chrom, index.Start, index.End, strings.Join(index.Transcripts, ","))
		}
	}
	log.Printf("Now you need run the command: 'samtools faidx %s'", outmRNA)
	defer reader.Close()
	defer writer.Close()
}

func NewPreGeneBasedCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "gb",
		Short: "Prepare Genebased database",
		Run: func(cmd *cobra.Command, args []string) {
			genome, _ := cmd.Flags().GetString("genome")
			refgene, _ := cmd.Flags().GetString("refgene")
			maneSelect, _ := cmd.Flags().GetString("mane_select")
			ncbiGeneInfo, _ := cmd.Flags().GetString("ncbi_gene_info")
			dbpath, _ := cmd.Flags().GetString("dbpath")
			builder, _ := cmd.Flags().GetString("builder")
			name, _ := cmd.Flags().GetString("name")
			step, _ := cmd.Flags().GetInt("step")
			if genome == "" || refgene == "" || dbpath == "" || builder == "" || name == "" || maneSelect == "" || ncbiGeneInfo == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				RunPreGeneBased(refgene, maneSelect, ncbiGeneInfo, genome, strings.ToLower(builder), step, path.Join(dbpath, builder, name))
			}
		},
	}
	cmd.Flags().StringP("genome", "g", "", "Reference Fasta File")
	cmd.Flags().StringP("refgene", "r", "", "RefGene File")
	cmd.Flags().StringP("dbpath", "d", "", "Database Directory")
	cmd.Flags().StringP("name", "n", "", "Database Name")
	cmd.Flags().StringP("builder", "b", "hg19", "Database Path")
	cmd.Flags().StringP("mane_select", "m", "", "Mane Select File, gzip")
	cmd.Flags().StringP("ncbi_gene_info", "c", "", "NCBI Gene Info file, gzip")
	cmd.Flags().IntP("step", "L", 300000, "Transcript Index Step Length")
	return cmd
}
