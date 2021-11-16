package del

import (
	"OpenAnno/anno/gene/snv/snp"
	"OpenAnno/db/transcript"
	"OpenAnno/variant"
	"strings"
)

type GeneAnnoItem struct {
	snp.GeneAnnoItem
}

func (a *GeneAnnoItem) AnnoInGene(del variant.Snv, refgene transcript.Transcript) {
	a.GeneSymbol = refgene.Gene
	a.GeneEntrezId = refgene.EntrezId
	a.Transcript = refgene.Transcript
	regionIndexes, lenL, lenR := refgene.Regions.FindMany(del.Start, del.End, refgene.Strand)
	regions := make(transcript.Regions, 0)
	for _, i := range regionIndexes {
		regions = append(regions, refgene.Regions[i])
	}
	if len(regionIndexes) == 1 {
		regionIndex := regionIndexes[0]
		region := refgene.Regions[regionIndex]
		if region.Type == "intron" {
			a.AnnoInIntron(del, refgene, regionIndex, lenL)
		} else if strings.HasPrefix(region.Type, "UTR") {
			a.AnnoInUTR(refgene, regionIndex)
		} else {
			a.AnnoInCDS(del, refgene, lenL, lenR, regions)
		}
	} else {
		if regions.HasCds() {
			a.AnnoInCDS(del, refgene, lenL, lenR, regions)
		} else {
			a.AnnoInMulti(regions)
		}
	}
}
