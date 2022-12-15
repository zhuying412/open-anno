package pre

import (
	"bufio"
	"fmt"
	"log"
	"open-anno/pkg"
	"os"
	"path"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type PreGeneParam struct {
	Gene2Refseq string `validate:"required,pathexists"`
	GeneInfo    string `validate:"required,pathexists"`
	GenePred    string `validate:"required,pathexists"`
	Output      string `validate:"required"`
}

func (this PreGeneParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this PreGeneParam) ReadGene2Refseq() (map[string]string, error) {
	transToId := make(map[string]string)
	reader, err := pkg.NewIOReader(this.Gene2Refseq)
	if err != nil {
		return transToId, err
	}
	defer reader.Close()
	scanner := pkg.NewCSVScanner(reader)
	for scanner.Scan() {
		row := scanner.Row()
		trans := strings.Split(row["RNA_nucleotide_accession.version"], ".")[0]
		entrezId := row["GeneID"]
		transToId[trans] = entrezId
	}
	return transToId, err
}

func (this PreGeneParam) ReadNCBIGeneInfo() (map[string]map[string]string, map[string]map[string]string, error) {
	symbolToId := make(map[string]map[string]string)
	synonymsToId := make(map[string]map[string]string)
	reader, err := pkg.NewIOReader(this.GeneInfo)
	if err != nil {
		return symbolToId, synonymsToId, err
	}
	defer reader.Close()
	scanner := pkg.NewCSVScanner(reader)
	for scanner.Scan() {
		row := scanner.Row()
		chrom := row["chromosome"]
		entrezId := row["GeneID"]
		symbol := row["Symbol"]
		synonyms := strings.Split(row["Synonyms"], "|")
		if _, ok := symbolToId[chrom]; !ok {
			symbolToId[chrom] = make(map[string]string)
			if _, ok := symbolToId[chrom][symbol]; !ok {
				symbolToId[chrom][symbol] = entrezId
			}
		}
		if _, ok := synonymsToId[chrom]; !ok {
			synonymsToId[chrom] = make(map[string]string)
			for _, name := range synonyms {
				if _, ok := synonymsToId[chrom][name]; !ok {
					synonymsToId[chrom][name] = entrezId
				}
			}
		}
	}
	return symbolToId, synonymsToId, err
}

func (this PreGeneParam) GetSymbolToId(transToId map[string]string, ncbiSymbolToId, ncbiSynonymsToId map[string]map[string]string) (map[string]map[string]string, error) {
	symbolToId := make(map[string]map[string]string)

	log.Printf("Read GenePred  from %s ...", this.GenePred)
	reader, err := pkg.NewIOReader(this.GenePred)
	if err != nil {
		return symbolToId, err
	}
	defer reader.Close()
	scanner := bufio.NewScanner(reader)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		chrom := fields[2]
		trans := strings.Split(fields[1], ".")[0]
		symbol := fields[12]
		if _, ok := symbolToId[chrom]; !ok {
			symbolToId[chrom] = make(map[string]string)
		}
		if _, ok := symbolToId[chrom][symbol]; !ok {
			if entrezId, ok := transToId[trans]; ok {
				symbolToId[chrom][symbol] = entrezId
			} else {
				if entrezId, ok = ncbiSymbolToId[chrom][symbol]; ok {
					symbolToId[chrom][symbol] = entrezId
				} else {
					if entrezId, ok = ncbiSynonymsToId[chrom][symbol]; ok {
						symbolToId[chrom][symbol] = entrezId
					}
				}
			}
		}
	}
	return symbolToId, err
}

func (this PreGeneParam) Run() error {
	log.Printf("Read NCBI Gene2Refseq from %s ...", this.Gene2Refseq)
	transToId, err := this.ReadGene2Refseq()
	if err != nil {
		return err
	}
	log.Printf("Read NCBI GeneInfo from %s ...", this.GeneInfo)
	ncbiSymbolToId, ncbiSynonymsToId, err := this.ReadNCBIGeneInfo()
	if err != nil {
		return err
	}
	symbolToId, err := this.GetSymbolToId(transToId, ncbiSymbolToId, ncbiSynonymsToId)
	if err != nil {
		return err
	}
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprint(writer, "Chrom\tSymbol\tEntrezId\n")
	for chrom, data := range symbolToId {
		for symbol, entrezId := range data {
			fmt.Fprintf(writer, "%s\t%s\t%s\n", chrom, symbol, entrezId)
		}
	}
	return nil
}

func NewGeneCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "gene",
		Short: "Get Gene Symbol to EntrezId",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreGeneParam
			param.Gene2Refseq, _ = cmd.Flags().GetString("gene2refseq")
			param.GeneInfo, _ = cmd.Flags().GetString("geneinfo")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.Output, _ = cmd.Flags().GetString("output")
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
	cmd.Flags().StringP("genepred", "p", "", "Input GenePred File")
	cmd.Flags().StringP("gene2refseq", "r", "", "Input NCBI Gene2Refseq File")
	cmd.Flags().StringP("geneinfo", "i", "", "Input NCBI GeneInfo File")
	cmd.Flags().StringP("output", "o", "", "Output File")
	return cmd
}
