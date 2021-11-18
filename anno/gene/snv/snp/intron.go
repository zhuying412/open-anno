package snp

import (
	"OpenAnno/pkg/transcript"
	"OpenAnno/pkg/variant"
	"fmt"
)

func (a *GeneAnnoItem) AnnoInIntron(snp variant.Snv, trans transcript.Transcript, regionIndex int, exonLen int) {
	region := trans.Regions[regionIndex]
	prevRegion, _ := trans.Regions.GetPrev(regionIndex, trans.Strand)
	nextRegion, _ := trans.Regions.GetNext(regionIndex, trans.Strand)
	distance1 := snp.Start - region.Start + 1
	distance2 := region.End - snp.End + 1
	//var closestRegion transcript.Region
	//var pos, distance int
	//var flag byte
	ref, alt := snp.Ref, snp.Alt
	distance := distance1
	closestRegion, pos, flag := prevRegion, exonLen, '+'
	if trans.Strand == '-' {
		ref.ReverseComplementing()
		alt.ReverseComplementing()
	}
	if distance1 > distance2 {
		distance = distance2
	}
	if (trans.Strand == '+') != (distance1 <= distance2) {
		closestRegion, pos, flag = nextRegion, exonLen+1, '-'
	}
	a.SetRegion("intronic")
	if distance <= SplicingDistance {
		a.SetRegion("splicing")
		if !closestRegion.IsCDS() {
			a.SetRegion(string(closestRegion.Type) + "_splicing")
		}
		a.SetEvent("splicing")
	}
	if closestRegion.IsCDS() {
		a.SetExon(closestRegion.ExonOrder)
		if snp.Type() == variant.SnvType_SNP {
			a.SetNAChange(fmt.Sprintf("c.%d%c%d%s>%s", pos, flag, distance, ref, alt))
		}
		if snp.Type() == variant.SnvType_INS {
			a.SetNAChange(fmt.Sprintf("c.%d%c%dins%s", pos, flag, distance, alt))
		}
	}
}
