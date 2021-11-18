package db

import (
	"OpenAnno/pkg/gene"
	"OpenAnno/pkg/utils"
	"github.com/spf13/viper"
	"log"
	"path"
	"strings"
)

var SymbolToEntrezId gene.EntrezIdMap

func InitGene() {
	if len(SymbolToEntrezId) == 0 {
		symbolToEntrezIdFile := path.Join(viper.GetString("db.path"), viper.GetString("db.gene_entrez"))
		log.Printf("read %s", symbolToEntrezIdFile)
		SymbolToEntrezId = readSymboToIdFile(symbolToEntrezIdFile)
	}
}

func readSymboToIdFile(symboToIdPath string) gene.EntrezIdMap {
	entrezIdMap := make(gene.EntrezIdMap)
	fi, reader := utils.OpenFile(symboToIdPath)
	defer utils.CloseFile(fi)
	for {
		line, isEof := utils.ReadLine(reader, 0)
		if isEof {
			break
		}
		fields := strings.Split(line, "\t")
		symbol, entrezId := fields[0], fields[1]
		entrezIdMap[symbol] = entrezId
	}
	return entrezIdMap
}
