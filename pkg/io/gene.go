package io

import "open-anno/pkg/scheme"

func ReadGeneInfo(infile string) (scheme.GeneInfoMap, error) {
	geneInfo := make(scheme.GeneInfoMap)
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
		gene := scheme.GeneInfo{
			Symbol:   row["Symbol"],
			EntrezId: row["EntrezId"],
			Chrom:    row["Chrom"],
		}
		if _, ok := geneInfo[gene.Chrom]; !ok {
			geneInfo[gene.Chrom] = make(map[string]scheme.GeneInfo)
		}
		geneInfo[gene.Chrom][gene.Symbol] = gene

	}
	return geneInfo, err
}
