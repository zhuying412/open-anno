package pre

import (
	"bufio"
	"errors"
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
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

func CheckPathExists(fl validator.FieldLevel) bool {
	path := fl.Field().String()
	_, err := os.Stat(path)
	return os.IsExist(err)
}

type PreGBParam struct {
	Genome       string `validate:"required,pathexists"`
	Refgene      string `validate:"required,pathexists"`
	Gene2Refseq  string `validate:"pathexists"`
	NcbiGeneInfo string `validate:"pathexists"`
	DBpath       string `validate:"required"`
	Builder      string `validate:"required"`
	Name         string `validate:"required"`
	IndexStep    int    `validate:"required"`
}

func (this PreGBParam) Outdir() string {
	outdir := path.Join(this.DBpath, this.Builder)
	if _, err := os.Stat(outdir); os.IsNotExist(err) {
		err := os.MkdirAll(outdir, os.ModePerm)
		if err != nil {
			log.Fatal(err)
		}
	}
	return outdir
}

func (this PreGBParam) OutmRNA() string {
	return path.Join(this.Outdir(), this.Name+"_mRNA.fa")
}

func (this PreGBParam) OutTransIndex() string {
	return path.Join(this.Outdir(), this.Name+".txt.idx")
}

func (this PreGBParam) OutRefGene() string {
	return path.Join(this.Outdir(), this.Name+".txt")
}

func (this PreGBParam) OutGene2Refseq() string {
	return path.Join(this.Outdir(), "Homo_sapiens.gene2refseq.gz")
}

func (this PreGBParam) OutNcbiGeneInfo() string {
	return path.Join(this.Outdir(), "Homo_sapiens.gene_info.gz")
}

func (this PreGBParam) NewGenomeFaidx() (*faidx.Faidx, error) {
	log.Printf("Read genome: %s ...", this.Genome)
	return faidx.New(this.Genome)
}

func (this PreGBParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	seq.SetGenome(this.Builder)
	return nil
}

func (this PreGBParam) WriteAndIndexmRNA(transcripts refgene.Transcripts) error {
	log.Printf("Write and Index mRNA: %s ...", this.OutmRNA())
	writer, err := io.NewIoWriter(this.OutmRNA())
	if err != err {
		return err
	}
	defer writer.Close()
	fai, err := this.NewGenomeFaidx()
	if err != err {
		return err
	}
	for _, trans := range transcripts {
		sequence, err := seq.Fetch(fai, trans.Chrom, trans.TxStart-1, trans.TxEnd)
		if err != nil {
			return err
		}
		sequence = strings.ToUpper(sequence)
		fmt.Fprintf(writer, ">%s:%s:%s\n%s\n", trans.Chrom, trans.Gene, trans.Name, sequence)
	}
	command := exec.Command("samtools", "faidx", this.OutmRNA())
	err = command.Run()
	if err != nil {
		log.Print(err)
		log.Printf("Now you need run the command: 'samtools faidx %s'", this.OutmRNA())
	}
	return nil
}

func (this PreGBParam) WriteTransIndex(transcripts refgene.Transcripts) error {
	log.Printf("Write Transcript Index: %s ...", this.OutTransIndex())
	writer, err := io.NewIoWriter(this.OutTransIndex())
	if err != nil {
		return err
	}
	defer writer.Close()
	transIndexes := refgene.NewTransIndexes(this.IndexStep)
	for _, index := range transIndexes {
		index.SetTranscripts(transcripts)
		if len(index.Transcripts) > 0 {
			fmt.Fprintf(writer, "%s\t%d\t%d\t%s\n", index.Chrom, index.Start, index.End, strings.Join(index.Transcripts, ","))
		}
	}
	return err
}

func (this PreGBParam) WriteGene2Refseq() error {
	if this.Gene2Refseq == "" {
		_, err := os.Stat(this.OutGene2Refseq())
		if os.IsNotExist(err) {
			return errors.New(fmt.Sprintf("Not Found: %s", this.OutGene2Refseq()))
		}
		return nil
	}
	log.Printf("Create NCBI Gene2Refseq to %s ...", this.OutGene2Refseq())
	reader, err := io.NewIoReader(this.Gene2Refseq)
	if err != nil {
		return err
	}
	defer reader.Close()
	writer, err := io.NewIoWriter(this.OutGene2Refseq())
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

func (this PreGBParam) WriteNcbiGeneInfo() error {
	if this.NcbiGeneInfo == "" {
		_, err := os.Stat(this.OutNcbiGeneInfo())
		if os.IsNotExist(err) {
			return errors.New(fmt.Sprintf("Not Found: %s", this.OutNcbiGeneInfo()))
		}
		return nil
	}
	log.Printf("Copy NCBI Gene Info to %s ...", this.OutNcbiGeneInfo())
	return io.CopyFile(this.NcbiGeneInfo, this.OutNcbiGeneInfo())
}

func (this PreGBParam) WriteAndReadRefgene() (refgene.Transcripts, error) {
	var transcripts refgene.Transcripts
	log.Printf("Create refgene to %s ...", this.OutRefGene())
	reader, err := io.NewIoReader(this.Refgene)
	if err != nil {
		return transcripts, err
	}
	defer reader.Close()
	writer, err := io.NewIoWriter(this.OutRefGene())
	if err != nil {
		return transcripts, err
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
	log.Printf("Read refgene: %s ...", this.OutRefGene())
	return refgene.ReadRefgene(this.OutRefGene())
}

func (this PreGBParam) Run() error {
	// gene2refseq
	err := this.WriteGene2Refseq()
	if err != nil {
		return err
	}
	// gene_info
	err = this.WriteNcbiGeneInfo()
	if err != nil {
		return err
	}
	// refgene
	transcripts, err := this.WriteAndReadRefgene()
	if err != nil {
		return err
	}
	// mRNA
	err = this.WriteAndIndexmRNA(transcripts)
	if err != nil {
		return err
	}
	// index
	err = this.WriteTransIndex(transcripts)
	if err != nil {
		return err
	}
	return nil
}

func NewPreGeneBasedCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "gene",
		Short: "Prepare Genebased database",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreGBParam
			param.Genome, _ = cmd.Flags().GetString("genome")
			param.Refgene, _ = cmd.Flags().GetString("refgene")
			param.Gene2Refseq, _ = cmd.Flags().GetString("gene2refseq")
			param.NcbiGeneInfo, _ = cmd.Flags().GetString("ncbi_gene_info")
			param.DBpath, _ = cmd.Flags().GetString("dbpath")
			param.Builder, _ = cmd.Flags().GetString("builder")
			param.Name, _ = cmd.Flags().GetString("name")
			param.IndexStep, _ = cmd.Flags().GetInt("step")
			err := param.Valid()
			if err != nil {
				cmd.Help()
				log.Fatal(err)
			}
			err = param.Run()
			if err != nil {
				log.Fatal(err)
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
