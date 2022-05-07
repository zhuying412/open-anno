package pre

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"open-anno/pkg/io/refgene"
	"open-anno/pkg/seq"
	"os"
	"os/exec"
	"path"
	"strings"

	"github.com/brentp/faidx"
	"github.com/spf13/cobra"
)

func writeGeneID(maneSelect, ncbiGeneInfo, refGene, outGeneId string) error {
	geneSymbolToId, err := io.NewGeneSymbolToId(maneSelect, ncbiGeneInfo, refGene)
	if err != nil {
		return err
	}
	writer, err := os.Create(outGeneId)
	if err != err {
		return err
	}
	writer.Close()
	for symbol, entrezId := range geneSymbolToId {
		fmt.Fprintf(writer, "%s\t%s\n", symbol, entrezId)
	}
	return err
}

func writemRNA(transcripts refgene.Transcripts, fai *faidx.Faidx, outmRNA string) error {
	writer, err := os.Create(outmRNA)
	if err != err {
		return err
	}
	writer.Close()
	for _, trans := range transcripts {
		sequence, err := seq.Fetch(fai, trans.Chrom, trans.TxStart-1, trans.TxEnd)
		if err != nil {
			return err
		}
		sequence = strings.ToUpper(sequence)
		fmt.Fprintf(writer, ">%s:%s:%s\n%s\n", trans.Chrom, trans.Gene, trans.Name, sequence)
	}
	return err
}

func writeTransIndex(transcripts refgene.Transcripts, indexStep int, outIndex string) error {
	writer, err := os.Create(outIndex)
	if err != nil {
		return err
	}
	defer writer.Close()
	transIndexes := refgene.NewTransIndexes(indexStep)
	for _, index := range transIndexes {
		index.SetTranscripts(transcripts)
		if len(index.Transcripts) > 0 {
			fmt.Fprintf(writer, "%s\t%d\t%d\t%s\n", index.Chrom, index.Start, index.End, strings.Join(index.Transcripts, ","))
		}
	}
	return err
}

func RunPreGeneBased(refGene string, maneSelect string, ncbiGeneInfo string, fasta string, builder string, indexStep int, outdir string) {
	// param
	log.Println("Init parameters ...")
	seq.SetGenome(builder)
	if _, err := os.Stat(outdir); os.IsNotExist(err) {
		err := os.MkdirAll(outdir, os.ModePerm)
		if err != nil {
			log.Fatal(err)
		}
	}
	// entrez_id
	outManeSelect := path.Join(outdir, "MANE.summary.txt.gz")
	log.Printf("Copy MANE Select to %s ...", outManeSelect)
	pkg.CopyFile(maneSelect, outManeSelect)
	outNcbiGeneInfo := path.Join(outdir, "Homo_sapiens.gene_info.gz")
	log.Printf("Copy NCBI Gene Info to %s ...", outNcbiGeneInfo)
	pkg.CopyFile(ncbiGeneInfo, outNcbiGeneInfo)
	// refgene
	outRefGene := path.Join(outdir, "refgene.txt")
	log.Printf("Copy refgene to %s ...", outRefGene)
	pkg.CopyFile(refGene, outRefGene)
	log.Printf("Read refgene: %s ...", outRefGene)
	transcripts, err := refgene.ReadRefgene(outRefGene)
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
	err = writemRNA(transcripts, fai, outmRNA)
	if err != nil {
		log.Fatal(err)
	}
	command := exec.Command("samtools", "faidx", outmRNA)
	err = command.Run()
	if err != nil {
		log.Print(err)
		log.Printf("Now you need run the command: 'samtools faidx %s'", outmRNA)
	}
	// index
	outIndex := path.Join(outdir, "refgene.idx")
	log.Printf("Write Transcript Index: %s ...", outIndex)
	err = writeTransIndex(transcripts, indexStep, outIndex)
	if err != nil {
		log.Fatal(err)
	}
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
