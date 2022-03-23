package genebased

import (
	"fmt"
	"open-anno/pkg/gene"
)

type SnvGeneBased struct {
	Gene       string `json:"gene"`
	GeneID     string `json:"gene_id"`
	Transcript string `json:"transcript"`
	Region     string `json:"region"`
	NAChange   string `json:"na_change"`
	AAChange   string `json:"aa_change"`
	Event      string `json:"event"`
	Region2    string `json:"region2"`
}

func NewSnvGeneBased(trans gene.Transcript, regions ...gene.Region) SnvGeneBased {
	anno := SnvGeneBased{
		Gene:       trans.Gene,
		GeneID:     trans.GeneID,
		Transcript: trans.Name,
		Event:      ".",
		NAChange:   ".",
		AAChange:   ".",
	}
	if len(regions) == 0 {
		anno.Region = "."
		anno.Region2 = "."
	} else {
		var region gene.Region
		if len(regions) == 1 {
			region = regions[0]
		} else {
			region1, region2 := regions[0], regions[1]
			if region1.Start > region2.Start {
				region1, region2 = region2, region1
			}
			if trans.Strand == "-" {
				region1, region2 = region2, region1
			}
			if region1.Exists() && region2.Exists() && region1.Name() != region2.Name() {
				anno.Region2 = fmt.Sprintf("%s_%s", region1.Name(), region2.Name())
				anno.Region = anno.Region2
			} else {
				region = region1
				if !region1.Exists() {
					region = region2
				}
			}
		}
		if region.Exists() {
			anno.Region2 = regions[0].Name()
			if regions[0].Type == gene.RType_UTR {
				anno.Region = regions[0].Name()
			} else if regions[0].Type == gene.RType_INTRON {
				anno.Region = "intronic"
			} else {
				anno.Region = "exonic"
			}
		}
	}
	return anno
}

type CnvGeneBased struct {
	Gene       string `json:"gene"`
	GeneID     string `json:"gene_id"`
	Transcript string `json:"transcript"`
	CDS        string `json:"cds"`
	Region     string `json:"utr3"`
	Strand     string `json:"strand"`
}

func NewCnvGeneBased(trans gene.Transcript) CnvGeneBased {
	anno := CnvGeneBased{
		Gene:       trans.Gene,
		GeneID:     trans.GeneID,
		Transcript: trans.Name,
		Strand:     trans.Strand,
		CDS:        ".",
		Region:     ".",
	}
	return anno
}
