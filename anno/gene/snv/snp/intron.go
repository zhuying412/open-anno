package snp

import (
	"OpenAnno/anno/gene/snv"
	"OpenAnno/db/transcript"
	snv2 "OpenAnno/variant"
	"fmt"
)

func (a *GeneAnnoItem) AnnoInIntron(snp snv2.Snv, trans transcript.Transcript, regionIndex int, exonLen int) {
	region := trans.Regions[regionIndex]
	prevRegion, _ := trans.Regions.GetPrev(regionIndex, trans.Strand)
	nextRegion, _ := trans.Regions.GetNext(regionIndex, trans.Strand)
	distance1 := snp.Start - region.Start + 1
	distance2 := region.End - snp.End + 1
	var closestRegion transcript.Region
	var pos, distance int
	var flag byte
	if distance1 <= distance2 {
		distance = distance1
	} else {
		distance = distance2
	}
	if (trans.Strand == '+') == (distance1 <= distance2) {
		closestRegion = prevRegion
		pos = exonLen
		flag = '+'
	} else {
		closestRegion = nextRegion
		pos = exonLen + 1
		flag = '-'
	}
	a.SetRegion("intronic")
	if distance <= snv.SplicingDistance {
		a.SetRegion("splicing")
		if !closestRegion.IsCDS() {
			a.SetRegion(closestRegion.Type + "_splicing")
		}
	}
	if closestRegion.IsCDS() {
		a.SetExon(closestRegion.ExonOrder)
		if snp.Type() == "snp" {
			a.SetNAChange(fmt.Sprintf("c.%d%c%d%s>%s", pos, flag, distance, snp.Ref, snp.Alt))
		} else {
			a.SetNAChange(fmt.Sprintf("c.%d%c%dins%s", pos, flag, distance, snp.Alt))
		}
	}
}
