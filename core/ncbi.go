package core

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

type NcbiGene struct {
	Symbol    map[string]int
	Synonyms  map[string]int
	SymbolFna map[string]int
}

func (ncbi NcbiGene) GetEntrezId(symbol string) int {
	if entrezId, ok := ncbi.Symbol[symbol]; ok {
		return entrezId
	}
	if entrezId, ok := ncbi.SymbolFna[symbol]; ok {
		return entrezId
	}
	if entrezId, ok := ncbi.Synonyms[symbol]; ok {
		return entrezId
	}
	return -1
}

func (ncbi *NcbiGene) Read(ncbiGeneInfoFile string) {
	ncbi.Symbol = make(map[string]int)
	ncbi.Synonyms = make(map[string]int)
	ncbi.SymbolFna = make(map[string]int)
	if fp, err := os.Open(ncbiGeneInfoFile); err == nil {
		defer fp.Close()
		var lines [][]byte
		if strings.HasSuffix(strings.ToLower(ncbiGeneInfoFile), ".gz") {
			if reader, _err := gzip.NewReader(fp); _err == nil {
				if content, err := ioutil.ReadAll(reader); err == nil {
					lines = bytes.Split(content, []byte{'\n'})
				}
			}
		} else {
			reader := bufio.NewReader(fp)
			if content, err := ioutil.ReadAll(reader); err == nil {
				lines = bytes.Split(content, []byte{'\n'})
			}
		}
		for _, line := range lines {
			line = bytes.TrimSpace(line)
			if len(line) == 0 || line[0] == '#' {
				continue
			}
			field := strings.Split(string(line), "\t")
			if entrezId, err := strconv.Atoi(field[1]); err == nil {
				symbol, synonyms, symbolFna := field[2], strings.Split(field[4], "|"), field[10]
				if symbol != "-" && symbol != "." {
					ncbi.Symbol[symbol] = entrezId
				}
				if symbolFna != "-" && symbolFna != "." {
					ncbi.SymbolFna[symbolFna] = entrezId
				}
				for _, synonym := range synonyms {
					if synonym != "-" && synonym != "." {
						ncbi.Synonyms[synonym] = entrezId
					}
				}
			}
		}
	} else {
		panic(err)
	}
}
