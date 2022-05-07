package io

import (
	"bufio"
	"compress/gzip"
	"open-anno/pkg"
	"open-anno/pkg/seq"
	"os"
	"strings"
)

func readManeSelect(infile string) (map[string]string, error) {
	transToId := make(map[string]string)
	fi, err := os.Open(infile)
	if err != nil {
		return transToId, err
	}
	defer fi.Close()
	reader, err := gzip.NewReader(fi)
	if err != nil {
		return transToId, err
	}
	defer reader.Close()
	scanner := NewCSVScanner(reader)
	for scanner.Scan() {
		row := scanner.Row()
		trans := strings.Split(row["RefSeq_nuc"], ".")[0]
		ensembl := strings.Split(row["Ensembl_nuc"], ".")[0]
		entrezId := strings.Split(row["#NCBI_GeneID"], ":")[1]
		transToId[trans] = entrezId
		transToId[ensembl] = entrezId
	}
	return transToId, err
}

func readNCBIGeneInfo(infile string) (map[string]string, map[string]string, error) {
	symbolToId := make(map[string]string)
	synonymsToId := make(map[string]string)
	fi, err := os.Open(infile)
	if err != nil {
		return symbolToId, synonymsToId, err
	}
	defer fi.Close()
	reader, err := gzip.NewReader(fi)
	if err != nil {
		return symbolToId, synonymsToId, err
	}
	defer reader.Close()
	scanner := NewCSVScanner(reader)
	for scanner.Scan() {
		row := scanner.Row()
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
	reader, err := os.Open(infile)
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

func NewGeneSymbolToId(maneSelect, ncbiGeneInfo, refgene string) (map[string]string, error) {
	symbolToId := make(map[string]string)
	maneTransToId, err := readManeSelect(maneSelect)
	if err != nil {
		return symbolToId, err
	}
	ncbiSymbolToId, ncbiSynonymsToId, err := readNCBIGeneInfo(ncbiGeneInfo)
	if err != nil {
		return symbolToId, err
	}
	reader, err := os.Open(refgene)
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
			if entrezId, ok := maneTransToId[trans]; ok {
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
