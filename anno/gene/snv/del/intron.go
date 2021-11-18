package del

import (
	"OpenAnno/anno/gene/snv/snp"
	"OpenAnno/pkg/transcript"
	"OpenAnno/pkg/variant"
	"fmt"
)

func (a *GeneAnnoItem) AnnoInIntron(del variant.Snv, trans transcript.Transcript, regionIndex int, exonLen int) {
	region := trans.Regions[regionIndex]
	prevRegion, _ := trans.Regions.GetPrev(regionIndex, trans.Strand)
	nextRegion, _ := trans.Regions.GetNext(regionIndex, trans.Strand)
	startDistance1, endDistance1 := del.Start-region.Start+1, del.End-region.Start+1
	startDistance2, endDistance2 := region.End-del.Start+1, region.End-del.End+1
	startDistance, endDistance := startDistance1, endDistance1
	if startDistance1 > endDistance2 {
		startDistance, endDistance = startDistance2, endDistance2
	}
	closestRegion, pos, flag := prevRegion, exonLen, '+'
	if (trans.Strand == '+') != (startDistance1 <= endDistance2) {
		closestRegion, pos, flag = nextRegion, exonLen+1, '-'
	}
	a.SetRegion("intronic")
	if startDistance <= snp.SplicingDistance || endDistance <= snp.SplicingDistance {
		a.SetRegion("splicing")
		if !closestRegion.IsCDS() {
			a.SetRegion(string(closestRegion.Type) + "_splicing")
		}
		a.SetEvent("splicing")
	}
	if closestRegion.IsCDS() {
		a.SetExon(closestRegion.ExonOrder)
		a.SetNAChange(fmt.Sprintf("c.%d%c%ddel", pos, startDistance, flag))
		if startDistance != endDistance {
			dis1, dis2 := startDistance, endDistance
			if (flag == '+' && dis1 > dis2) || (flag == '-' && dis1 < dis2) {
				dis1, dis2 = endDistance, startDistance
			}
			a.SetNAChange(fmt.Sprintf("c.%d%c%d_%d%c%ddel", pos, flag, dis1, flag, dis2))
		}
	}
}
