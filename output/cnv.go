package output

import (
	"grandanno/gene_based/cnv"
	"grandanno/input"
)

type CnvOutput struct {
	Cnv       input.Cnv       `json:"Cnv"`
	GeneAnnos []cnv.GeneAnno  `json:"gene_anno"`
	OtherInfo input.OtherInfo `json:"other_info"`
}

type CnvOutputs []CnvOutput

func (c CnvOutputs) JsonLines() (lines []string) {
	for _, output := range c {
		lines = append(lines, ConvertToJSON(output))
	}
	return lines
}
