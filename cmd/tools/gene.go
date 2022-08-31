package tools

import (
	"bufio"
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"open-anno/pkg/seq"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type GeneInfoParam struct {
	Gene2Refseq string `validate:"required"`
	GeneInfo    string `validate:"required"`
	GenePred    string `validate:"required"`
	Ouput       string `validate:"required"`
}

func (this GeneInfoParam) Valid() error {
	validate := validator.New()
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return nil
}

func (this GeneInfoParam) ReadGene2Refseq() (map[string]string, error) {
	transToId := make(map[string]string)
	reader, err := io.NewIoReader(this.Gene2Refseq)
	if err != nil {
		return transToId, err
	}
	defer reader.Close()
	scanner := io.NewCSVScanner(reader)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return transToId, err
		}
		trans := strings.Split(row["RNA_nucleotide_accession.version"], ".")[0]
		entrezId := row["GeneID"]
		transToId[trans] = entrezId
	}
	return transToId, err
}

func (this GeneInfoParam) ReadNCBIGeneInfo() (map[string]map[string]string, map[string]map[string]string, error) {
	symbolToId := make(map[string]map[string]string)
	synonymsToId := make(map[string]map[string]string)
	reader, err := io.NewIoReader(this.GeneInfo)
	if err != nil {
		return symbolToId, synonymsToId, err
	}
	defer reader.Close()
	scanner := io.NewCSVScanner(reader)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return symbolToId, synonymsToId, err
		}
		chrom := pkg.FormatChrom(row["chromosome"])
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

func (this GeneInfoParam) GetSymbolToId() (map[string]map[string]string, error) {
	symbolToId := make(map[string]map[string]string)
	transToId, err := this.ReadGene2Refseq()
	if err != nil {
		return symbolToId, err
	}
	ncbiSymbolToId, ncbiSynonymsToId, err := this.ReadNCBIGeneInfo()
	if err != nil {
		return symbolToId, err
	}
	reader, err := io.NewIoReader(this.GenePred)
	if err != nil {
		return symbolToId, err
	}
	defer reader.Close()
	scanner := bufio.NewScanner(reader)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		chrom := pkg.FormatChrom(fields[2])
		if _, ok := seq.GENOME[chrom]; ok {
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
	}
	return symbolToId, err
}

func (this GeneInfoParam) Run() error {
	symbolToId, err := this.GetSymbolToId()
	if err != nil {
		return err
	}
	writer, err := io.NewIoWriter(this.Ouput)
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

func NewGeneInfoCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "geneinfo",
		Short: "Get Gene Symbol to EntrezId",
		Run: func(cmd *cobra.Command, args []string) {
			var param GeneInfoParam
			param.Gene2Refseq, _ = cmd.Flags().GetString("gene2refseq")
			param.GeneInfo, _ = cmd.Flags().GetString("geneinfo")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
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
	cmd.Flags().StringP("genepred", "p", "", "GenePred File")
	cmd.Flags().StringP("gene2refseq", "r", "", "Gene2Refseq File")
	cmd.Flags().StringP("geneinfo", "i", "", "GeneInfo File")
	return cmd
}
