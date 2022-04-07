package genebased

import (
	"fmt"
	"open-anno/pkg/gene"
	"sort"
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
	}
	nregions := make(gene.Regions, 0)
	for _, region := range regions {
		if region.Exists() {
			nregions = append(nregions, region)
		}
	}
	if trans.Strand == "+" {
		sort.Sort(nregions)
	} else {
		sort.Sort(sort.Reverse(nregions))
	}
	if len(nregions) > 0 {
		region1, region2 := nregions[0], nregions[len(nregions)-1]
		anno.Region2 = region1.Name()
		if region1.Name() != region2.Name() {
			anno.Region2 = fmt.Sprintf("%s_%s", region1.Name(), region2.Name())
		}
		if region1.Type == gene.RType_CDS || region2.Type == gene.RType_CDS {
			anno.Region = "exonic"
		} else {
			if region1.Type == gene.RType_UTR {
				anno.Region = region1.Name()
			} else {
				if region2.Type == gene.RType_UTR {
					anno.Region = region2.Name()
				} else {
					anno.Region = "intronic"
				}
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
	Region     string `json:"region"`
	Strand     string `json:"strand"`
	Position   string `json:"position"`
}

func NewCnvGeneBased(trans gene.Transcript) CnvGeneBased {
	anno := CnvGeneBased{
		Gene:       trans.Gene,
		GeneID:     trans.GeneID,
		Transcript: trans.Name,
		Strand:     trans.Strand,
		CDS:        ".",
		Region:     ".",
		Position:   fmt.Sprintf("%d-%d", trans.TxStart, trans.TxEnd),
	}
	if anno.GeneID == "" {
		anno.GeneID = "."
	}
	return anno
}
