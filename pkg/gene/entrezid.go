package gene

import (
	"bufio"
	"compress/gzip"
	"open-anno/pkg"
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
	scanner := bufio.NewScanner(reader)
	scanner.Scan()
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		trans := strings.Split(fields[5], ".")[0]
		ensembl := strings.Split(fields[7], ".")[0]
		entrezId := strings.Split(fields[0], ":")[1]
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
	scanner := bufio.NewScanner(reader)
	scanner.Scan()
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		symbol := fields[1]
		entrezId := fields[2]
		synonyms := strings.Split(fields[4], "|")
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

func NewGeneSymbolToId(maneSelect string, ncbiGeneInfo string, refgene string) (map[string]string, error) {
	symbolToId := make(map[string]string)
	maneTransToId, err := readManeSelect(maneSelect)
	if err != nil {
		return symbolToId, err
	}
	ncbiSymbolToId, ncbiSynonymsToId, err := readNCBIGeneInfo(ncbiGeneInfo)
	if err != nil {
		return symbolToId, err
	}
	refgeneTransToSymbols, err := readRefgene(refgene)
	if err != nil {
		return symbolToId, err
	}
	for trans, symbols := range refgeneTransToSymbols {
		transShort := strings.Split(trans, ".")[0]
		for _, symbol := range symbols {
			if _, ok := symbolToId[symbol]; ok {
				continue
			}
			if entrezId, ok := maneTransToId[transShort]; ok {
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

func ReadGeneSymbolToId(infile string) (map[string]string, error) {
	symbolToId := make(map[string]string)
	reader, err := os.Open(infile)
	if err != nil {
		return symbolToId, err
	}
	defer reader.Close()
	scanner := bufio.NewScanner(reader)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		symbol := fields[0]
		entrezId := fields[1]
		symbolToId[symbol] = entrezId
	}
	return symbolToId, err
}
