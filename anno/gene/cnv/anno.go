package cnv

import (
	"OpenAnno/anno"
	"OpenAnno/db/transcript"
)

type GeneAnnoItem struct {
	GeneSymbol   string `json:"gene_symbol"`
	GeneEntrezId string `json:"gene_entrez_id"`
	Transcript   string `json:"transcript"`
	Region       string `json:"region"`
	StartExon    int    `json:"StartExon"`
	EndExon      int    `json:"end_exon"`
}

func (a *GeneAnnoItem) SetGene(trans transcript.Transcript) {
	a.GeneSymbol = trans.Gene
	a.GeneEntrezId = trans.EntrezId
	a.Transcript = trans.Transcript
}

func (a *GeneAnnoItem) SetExon(exon int) {
	if exon < a.StartExon {
		a.StartExon = exon
	}
	if exon > a.EndExon {
		a.EndExon = exon
	}
}

func (a *GeneAnnoItem) SetRegion(region string) {
	a.Region = region
}

func (a GeneAnnoItem) InCodingOrSplicingRegion() bool {
	return a.StartExon > 0 && a.EndExon > 0
}

type GeneAnno []GeneAnnoItem

func (a GeneAnno) AnnoType() anno.AnnoType {
	return anno.AnnoType_GENE
}

func (a GeneAnno) HasCodingorSplicingRegion() bool {
	for _, _anno := range a {
		if _anno.InCodingOrSplicingRegion() {
			return true
		}
	}
	return false
}
