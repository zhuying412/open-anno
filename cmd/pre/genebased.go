package pre

import (
	"bufio"
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"open-anno/pkg/scheme"
	"open-anno/pkg/seq"
	"os"
	"path"
	"strings"

	"github.com/brentp/faidx"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

func CheckPathExists(fl validator.FieldLevel) bool {
	path := fl.Field().String()
	_, err := os.Stat(path)
	return !os.IsNotExist(err)
}

type PreGBParam struct {
	Genome    string `validate:"required,pathexists"`
	GenePred  string `validate:"required,pathexists"`
	DBpath    string `validate:"required"`
	Builder   string `validate:"required"`
	Name      string `validate:"required"`
	IndexStep int    `validate:"required"`
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

func (this PreGBParam) OutGenePred() string {
	return path.Join(this.Outdir(), this.Name+".geneinfo.txt")
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

func (this PreGBParam) CreateAndIndexmRNA(transcripts scheme.Transcripts) error {
	log.Printf("Write and Index mRNA: %s ...", this.OutmRNA())
	log.Printf("Read genome: %s ...", this.Genome)
	fai, err := faidx.New(this.Genome)
	if err != err {
		return err
	}
	return io.CreateAndIndexmRNA(transcripts, fai, this.OutmRNA())
}

func (this PreGBParam) CreateTransIndex(transcripts scheme.Transcripts) error {
	log.Printf("Write Transcript Index: %s ...", this.OutTransIndex())
	return io.CreateTransIndexes(transcripts, this.IndexStep, this.OutTransIndex())
}

func (this PreGBParam) CreateAndReadRefgene() (scheme.Transcripts, error) {
	var transcripts scheme.Transcripts
	log.Printf("Create refgene to %s ...", this.OutGenePred())
	reader, err := io.NewIoReader(this.GenePred)
	if err != nil {
		return transcripts, err
	}
	defer reader.Close()
	writer, err := io.NewIoWriter(this.OutGenePred())
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
	log.Printf("Read refgene: %s ...", this.OutGenePred())
	return io.ReadGenePred(this.OutGenePred())
}

func (this PreGBParam) Run() error {
	// refgene
	transcripts, err := this.CreateAndReadRefgene()
	if err != nil {
		return err
	}
	// mRNA
	err = this.CreateAndIndexmRNA(transcripts)
	if err != nil {
		return err
	}
	// index
	err = this.CreateTransIndex(transcripts)
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
			param.GenePred, _ = cmd.Flags().GetString("genepred")
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
	cmd.Flags().StringP("genepred", "r", "", "RefGene File")
	cmd.Flags().StringP("dbpath", "d", "", "Database Directory")
	cmd.Flags().StringP("name", "n", "", "Database Name")
	cmd.Flags().StringP("builder", "b", "hg19", "Database Path")
	cmd.Flags().IntP("step", "L", 300000, "Transcript Index Step Length")
	return cmd
}
