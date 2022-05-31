package pre

import (
	"bufio"
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
	writer, err := io.NewIoWriter(outGeneId)
	if err != err {
		return err
	}
	defer writer.Close()
	for symbol, entrezId := range geneSymbolToId {
		fmt.Fprintf(writer, "%s\t%s\n", symbol, entrezId)
	}
	return err
}

func writemRNA(transcripts refgene.Transcripts, fai *faidx.Faidx, outmRNA string) error {
	writer, err := io.NewIoWriter(outmRNA)
	if err != err {
		return err
	}
	defer writer.Close()
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
	writer, err := io.NewIoWriter(outIndex)
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

func writeGene2Refseq(infile, outfile string) error {
	reader, err := io.NewIoReader(infile)
	if err != nil {
		return err
	}
	defer reader.Close()
	writer, err := io.NewIoWriter(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	scanner := bufio.NewScanner(reader)
	for scanner.Scan() {
		text := scanner.Text()
		if strings.HasPrefix(text, "#tax_id") || strings.HasPrefix(text, "9606") {
			fmt.Fprintf(writer, "%s\n", text)
		}
	}
	return err
}

func writeRefgene(infile, outfile string) error {
	reader, err := io.NewIoReader(infile)
	if err != nil {
		return err
	}
	defer reader.Close()
	writer, err := io.NewIoWriter(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	scanner := bufio.NewScanner(reader)
	for scanner.Scan() {
		text := scanner.Text()
		chrom := pkg.FormatChrom(strings.Split(text, "\t")[1])
		if _, ok := seq.GENOME[chrom]; ok {
			fmt.Fprintf(writer, "%s\n", text)
		}
		chrom = pkg.FormatChrom(strings.Split(text, "\t")[2])
		if _, ok := seq.GENOME[chrom]; ok {
			fmt.Fprintf(writer, "%s\n", text)
		}
	}
	return err
}

func RunPreGeneBased(dbpath, name, refGene, gene2refseq, ncbiGeneInfo, fasta, builder string, indexStep int) {
	// param
	log.Println("Init parameters ...")
	seq.SetGenome(builder)
	builderDir := path.Join(dbpath, builder)
	outdir := path.Join(builderDir, name)
	if _, err := os.Stat(outdir); os.IsNotExist(err) {
		err := os.MkdirAll(outdir, os.ModePerm)
		if err != nil {
			log.Fatal(err)
		}
	}
	if gene2refseq != "" {
		// gene2refseq
		outGene2Refseq := path.Join(builderDir, "Homo_sapiens.gene2refseq.gz")
		log.Printf("Create NCBI Gene2Refseq to %s ...", outGene2Refseq)
		writeGene2Refseq(gene2refseq, outGene2Refseq)
	}
	if ncbiGeneInfo != "" {
		// gene_info
		outNcbiGeneInfo := path.Join(builderDir, "Homo_sapiens.gene_info.gz")
		log.Printf("Copy NCBI Gene Info to %s ...", outNcbiGeneInfo)
		io.CopyFile(ncbiGeneInfo, outNcbiGeneInfo)
	}
	// refgene
	outRefGene := path.Join(outdir, "GenePred.txt")
	log.Printf("Create refgene to %s ...", outRefGene)
	err := writeRefgene(refGene, outRefGene)
	if err != nil {
		log.Fatal(err)
	}
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
	outIndex := path.Join(outdir, "GenePred.idx")
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
			gene2refseq, _ := cmd.Flags().GetString("gene2refseq")
			ncbiGeneInfo, _ := cmd.Flags().GetString("ncbi_gene_info")
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
				RunPreGeneBased(dbpath, name, refgene, gene2refseq, ncbiGeneInfo, genome, strings.ToLower(builder), step)
			}
		},
	}
	cmd.Flags().StringP("genome", "g", "", "Reference Fasta File")
	cmd.Flags().StringP("refgene", "r", "", "RefGene File")
	cmd.Flags().StringP("dbpath", "d", "", "Database Directory")
	cmd.Flags().StringP("name", "n", "", "Database Name")
	cmd.Flags().StringP("builder", "b", "hg19", "Database Path")
	cmd.Flags().StringP("gene2refseq", "m", "", "NCBI Gene to Refseq file, gzip")
	cmd.Flags().StringP("ncbi_gene_info", "c", "", "NCBI Gene Info file, gzip")
	cmd.Flags().IntP("step", "L", 300000, "Transcript Index Step Length")
	return cmd
}
