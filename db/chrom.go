package db

import (
	"OpenAnno/pkg/gene"
	"OpenAnno/pkg/utils"
	"github.com/spf13/viper"
	"log"
	"path"
	"strconv"
	"strings"
)

var ChromList gene.Chromosomes

func InitChrom() {
	if len(ChromList) == 0 {
		refDictFile := path.Join(viper.GetString("db.path"), viper.GetString("db.reference_dict"))
		log.Printf("read %s", refDictFile)
		ChromList = readReferenceDictFile(refDictFile)
	}
}

func readReferenceDictFile(referenceDictPath string) gene.Chromosomes {
	chroms := make(gene.Chromosomes, 0)
	fi, reader := utils.OpenFile(referenceDictPath)
	defer utils.CloseFile(fi)
	for {
		line, isEof := utils.ReadLine(reader, 0)
		if isEof {
			break
		}
		fields := strings.Split(line, "\t")
		if fields[0] == "@SQ" {
			name := strings.Split(string(fields[1]), ":")[1]
			length, err := strconv.Atoi(strings.Split(string(fields[2]), ":")[1])
			if err != nil {
				log.Panic(err)
			}
			chroms = append(chroms, gene.Chromosome{Name: name, Length: length})
		}
	}
	return chroms
}
