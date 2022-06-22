package io

import (
	"bufio"
	"open-anno/pkg"
	"open-anno/pkg/seq"
	"strings"
)

func readGene2Refseq(infile string) (map[string]string, error) {
	transToId := make(map[string]string)
	reader, err := NewIoReader(infile)
	if err != nil {
		return transToId, err
	}
	defer reader.Close()
	scanner := NewCSVScanner(reader)
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

func readNCBIGeneInfo(infile string) (map[string]string, map[string]string, error) {
	symbolToId := make(map[string]string)
	synonymsToId := make(map[string]string)
	reader, err := NewIoReader(infile)
	if err != nil {
		return symbolToId, synonymsToId, err
	}
	defer reader.Close()
	scanner := NewCSVScanner(reader)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return symbolToId, synonymsToId, err
		}
		entrezId := row["GeneID"]
		symbol := row["Symbol"]
		synonyms := strings.Split(row["Synonyms"], "|")
		if _, ok := symbolToId[symbol]; !ok {
			symbolToId[symbol] = entrezId
		}
		for _, name := range synonyms {
			if _, ok := synonymsToId[name]; !ok {
				synonymsToId[name] = entrezId
			}
		}
	}
	return symbolToId, synonymsToId, err

}

func readRefgene(infile string) (map[string][]string, error) {
	transToSymbols := make(map[string][]string)
	reader, err := NewIoReader(infile)
	if err != nil {
		return transToSymbols, err
	}
	defer reader.Close()
	scanner := bufio.NewScanner(reader)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		trans := fields[1]
		symbol := fields[12]
		if genes, ok := transToSymbols[trans]; ok {
			if pkg.FindArr(genes, symbol) < 0 {
				genes = append(genes, symbol)
				transToSymbols[trans] = genes
			}
		} else {
			transToSymbols[trans] = []string{symbol}
		}
	}
	return transToSymbols, err
}

func NewGeneSymbolToId(gene2refseq, ncbiGeneInfo, refgene string) (map[string]string, error) {
	symbolToId := make(map[string]string)
	transToId, err := readGene2Refseq(gene2refseq)
	if err != nil {
		return symbolToId, err
	}
	ncbiSymbolToId, ncbiSynonymsToId, err := readNCBIGeneInfo(ncbiGeneInfo)
	if err != nil {
		return symbolToId, err
	}
	reader, err := NewIoReader(refgene)
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
			if entrezId, ok := transToId[trans]; ok {
				symbolToId[symbol] = entrezId
			} else {
				if entrezId, ok = ncbiSymbolToId[symbol]; ok {
					symbolToId[symbol] = entrezId
				} else {
					if entrezId, ok = ncbiSynonymsToId[symbol]; ok {
						symbolToId[symbol] = entrezId
					}
				}
			}
		}
	}
	return symbolToId, err
}
