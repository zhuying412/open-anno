package pre

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"os"
	"path"
	"strconv"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type PreGenePredParam struct {
	Input       string `validate:"required,pathexists"`
	GenomeIndex string `validate:"required,pathexists"`
	MaxLength   int    `validate:"required"`
	Output      string `validate:"required"`
	IndexStep   int    `validate:"required"`
}

func (this PreGenePredParam) OutputIndex() string {
	return this.Output + ".idx"
}

func (this PreGenePredParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

// ReadGenomeIndex 读取Genome faidx文件，如hg38.fasta.fai
func (this PreGenePredParam) ReadGenomeIndex() (map[string]int, error) {
	chromLen := make(map[string]int)
	reader, err := os.Open(this.GenomeIndex)
	if err != nil {
		return chromLen, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		text := scanner.Text()
		fields := strings.Split(text, "\t")
		chrom := fields[0]
		if len(chrom) > this.MaxLength {
			continue
		}
		length, err := strconv.Atoi(strings.Split(fields[2], ":")[1])
		if err != nil {
			return chromLen, err
		}
		chromLen[chrom] = length
	}
	return chromLen, nil
}

func (this PreGenePredParam) CreateGenePred(chromLen map[string]int) error {
	reader, err := pkg.NewIOReader(this.Input)
	if err != nil {
		return err
	}
	defer reader.Close()
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		if _, ok := chromLen[row[2]]; ok {
			fmt.Fprintln(writer, scanner.Text())
		}
	}
	return nil
}

func (this PreGenePredParam) Run() error {
	// 读取Genome Index
	log.Printf("Read Genome Index: %s ...", this.GenomeIndex)
	chromLen, err := this.ReadGenomeIndex()
	if err != err {
		return err
	}
	// 创建 GenePred
	if this.Input != this.Output {
		log.Printf("Create GenePred: %s ...", this.Output)
		this.CreateGenePred(chromLen)
	}
	// 构建GenePred索引
	log.Printf("Read GenePred: %s ...", this.Input)
	gpes, err := pkg.ReadGenePred(this.Input)
	if err != err {
		return err
	}
	log.Printf("Create TransIndex: %s ...", this.OutputIndex())
	return pkg.CreateTransIndexes(gpes, chromLen, this.IndexStep, this.OutputIndex())
}

func NewPreGenePredCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "gpe",
		Short: "Prepare GenePred Files",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreGenePredParam
			param.GenomeIndex, _ = cmd.Flags().GetString("genome_index")
			param.MaxLength, _ = cmd.Flags().GetInt("max_len")
			param.Input, _ = cmd.Flags().GetString("input")
			param.Output, _ = cmd.Flags().GetString("output")
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
	cmd.Flags().StringP("genome_index", "g", "", "Input Reference Faidx File")
	cmd.Flags().StringP("input", "i", "", "Input GenePred File")
	cmd.Flags().StringP("output", "o", "", "Output GenePred File")
	cmd.Flags().IntP("max_len", "m", 5, "Max Length of Chromosome Name")
	cmd.Flags().IntP("step", "l", 300000, "Transcript Index Step Length")
	return cmd
}
