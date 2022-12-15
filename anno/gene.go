package anno

import (
	"open-anno/pkg"
	"os"
	"strconv"
	"strings"
)

func ReadGene(infile string) (map[string]map[string]string, error) {
	gene := make(map[string]map[string]string)
	reader, err := pkg.NewIOReader(infile)
	if err != nil {
		return gene, err
	}
	defer reader.Close()
	scanner := pkg.NewCSVScanner(reader)
	scanner.Scan()
	for scanner.Scan() {
		row := scanner.Row()
		chrom := row["Chrom"]
		entrezId := row["EntrezId"]
		symbol := row["Symbol"]
		if _, ok := gene[chrom]; !ok {
			gene[chrom] = make(map[string]string)
		}
		gene[chrom][symbol] = entrezId

	}
	return gene, nil
}

// ReadGenomeDict， 读取GenomeDict文件，如hg38.dict
func ReadGenomeDict(refDictFile string) (map[string]int, error) {
	genome := make(map[string]int)
	reader, err := os.Open(refDictFile)
	if err != nil {
		return genome, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		text := scanner.Text()
		fields := strings.Split(text, "\t")
		if fields[0] == "@SQ" {
			chrom := strings.Split(fields[1], ":")[1]
			length, err := strconv.Atoi(strings.Split(fields[2], ":")[1])
			if err != nil {
				return genome, err
			}
			genome[chrom] = length
		}
	}
	return genome, nil
}
