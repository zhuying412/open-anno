package anno

import (
	"open-anno/pkg"
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

type AnnoInfo struct {
	Key   string `json:"key"`
	Value string `json:"value"`
}

type AnnoInfos map[string][]AnnoInfo

func (this *AnnoInfos) Add(pk, key, val string) {
	if items, ok := (*this)[pk]; ok {
		(*this)[pk] = append(items, AnnoInfo{Key: key, Value: val})
	} else {
		(*this)[pk] = []AnnoInfo{{Key: key, Value: val}}
	}
}
