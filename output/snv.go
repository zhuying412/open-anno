package output

import (
	"grandanno/gene_based/snv"
	"grandanno/input"
)

type SnvOutput struct {
	Snv       input.Snv       `json:"snv"`
	GeneAnnos []snv.GeneAnno  `json:"gene_anno"`
	OtherInfo input.OtherInfo `json:"other_info"`
}

type SnvOutputs []SnvOutput

func (s SnvOutputs) JsonLines() (lines []string) {
	for _, output := range s {
		lines = append(lines, ConvertToJSON(output))
	}
	return lines
}
