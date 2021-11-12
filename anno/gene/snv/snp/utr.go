package snp

import (
	"OpenAnno/db/transcript"
)

func (a *GeneAnnoItem) AnnoInUTR(trans transcript.Transcript, regionIndex int) {
	region := trans.Regions[regionIndex]
	a.SetRegion(region.Type)
}
