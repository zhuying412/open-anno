package snv

import (
	"OpenAnno/anno"
	"OpenAnno/db/transcript"
	"OpenAnno/variant"
)

type IGeneAnnoItem interface {
	InCodingOrSplicingRegion() bool
	SetExon(exon int)
	SetGene(trans transcript.Transcript)
	SetRegion(region string)
	SetEvent(event string)
	SetNAChange(change string)
	SetAAChange(change string)
	AnnoInGene(snv variant.Snv, transcript transcript.Transcript)
}

type GeneAnno []IGeneAnnoItem

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
