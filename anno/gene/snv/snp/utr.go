package snp

import "OpenAnno/pkg/transcript"

func (a *GeneAnnoItem) AnnoInUTR(trans transcript.Transcript, regionIndex int) {
	region := trans.Regions[regionIndex]
	a.SetRegion(string(region.Type))
}
