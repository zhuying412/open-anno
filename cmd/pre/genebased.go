package pre

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/gene"
	"os"
	"path"

	"github.com/brentp/faidx"
	"github.com/spf13/cobra"
)

func PreGeneBased(refgene string, fasta string, builder string, indexStep int, outdir string) {
	log.Println("Init parameters ...")
	genome := gene.GENOME_HG19
	if builder == "hg38" {
		genome = gene.GENOME_HG19
	}
	if _, err := os.Stat(outdir); os.IsNotExist(err) {
		err := os.MkdirAll(outdir, os.ModePerm)
		if err != nil {
			log.Fatal(err)
		}
	}
	log.Printf("Read refgene: %s ...", refgene)
	transcripts, err := gene.ReadRefgene(refgene)
	if err != nil {
		log.Fatal(err)
	}
	log.Printf("Read genome: %s ...", fasta)
	fai, err := faidx.New(fasta)
	if err != nil {
		log.Fatal(err)
	}
	log.Printf("Write transcript: %s ...", outdir)
	writers := make(map[string]*os.File)
	for _, trans := range transcripts {
		err = trans.SetRegions(fai)
		if err != nil {
			log.Fatal(err)
		}
		if _, ok := genome[trans.Chrom]; ok {
			if _, ok := writers[trans.Chrom]; !ok {
				outfile := path.Join(outdir, fmt.Sprintf("chr%s.json", trans.Chrom))
				writers[trans.Chrom], err = os.Create(outfile)
				if err != nil {
					log.Fatal(err)
				}
			}
			writers[trans.Chrom].WriteString(pkg.ToJSON(trans) + "\n")
		}
	}
	for _, writer := range writers {
		err := writer.Close()
		if err != nil {
			log.Panic()
		}
	}
	log.Printf("Write transcript index: %s ...", outdir)
	writers = make(map[string]*os.File)
	transIndexes := gene.NewTransIndexes(genome, indexStep)
	for _, idx := range transIndexes {
		idx.SetTranscripts(transcripts)
		if _, ok := genome[idx.Chrom]; ok {
			if _, ok := writers[idx.Chrom]; !ok {
				outfile := path.Join(outdir, fmt.Sprintf("chr%s.idx.json", idx.Chrom))
				writers[idx.Chrom], err = os.Create(outfile)
				if err != nil {
					log.Fatal(err)
				}
			}
			writers[idx.Chrom].WriteString(pkg.ToJSON(idx) + "\n")
		}
	}
	for _, writer := range writers {
		err := writer.Close()
		if err != nil {
			log.Panic()
		}
	}
}

func NewPreGeneBasedCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "GeneBased",
		Short: "Prepare required transcript files",
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
				// fmt.Println(refgene, genome, builder, step, path.Join(dbpath, builder, name), dbpath)
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
