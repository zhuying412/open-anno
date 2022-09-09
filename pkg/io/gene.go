package io

import "open-anno/pkg/schema"

func ReadGeneInfo(infile string) (schema.GeneInfoMap, error) {
	geneInfo := make(schema.GeneInfoMap)
	reader, err := NewIoReader(infile)
	if err != nil {
		return geneInfo, err
	}
	defer reader.Close()
	scanner := NewCSVScanner(reader)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return geneInfo, err
		}
		gene := schema.GeneInfo{
			Symbol:   row["Symbol"],
			EntrezId: row["EntrezId"],
			Chrom:    row["Chrom"],
		}
		if _, ok := geneInfo[gene.Chrom]; !ok {
			geneInfo[gene.Chrom] = make(map[string]schema.GeneInfo)
		}
		geneInfo[gene.Chrom][gene.Symbol] = gene

	}
	return geneInfo, err
}
