package snp

import (
	"OpenAnno/db/transcript"
	"OpenAnno/seq"
	"OpenAnno/variant"
	"fmt"
)

func (a *GeneAnnoItem) AnnoInCDS(snp variant.Snv, trans transcript.Transcript, regionIndex int, exonLen int) {
	region := trans.Regions[regionIndex]
	var pos int
	a.SetExon(region.ExonOrder)
	a.SetRegion("exonic")
	if trans.Strand == '+' {
		pos = exonLen + snp.Start - region.Start + 1
	} else {
		pos = exonLen + region.End - snp.Start + 1
	}
	if trans.IsCmpl() {
		cdna, protein := trans.Cdna, trans.Protein
		newCdna := cdna.ChangeWithSnp(pos, snp.Alt.Base(0))
		newProtein := newCdna.Translate(snp.Chrom == "MT")
		for i := 0; i < cdna.Len(); i++ {
			if j, na1, na2 := i/3, cdna.Base(i), newCdna.Base(i); na1 != na2 && j < protein.Len() {
				aa1, aa2 := protein.Base(j), newProtein.Base(j)
				if aa1 == aa2 {
					a.SetEvent("snp_synonymous")
				} else {
					if aa1 == '*' {
						a.SetEvent("snp_stoploss")
					} else if aa2 == '*' {
						a.SetEvent("snp_stopgain")
					} else {
						if j == 0 {
							a.SetEvent("snp_startloss")
						}
						a.SetEvent("snp_nonsynonymous")
					}
				}
				a.SetNAChange(fmt.Sprintf("c.%d%c>%c", i+1, na1, na2))
				a.SetAAChange(fmt.Sprintf("p.%s%d%s", seq.AAMap[aa1], j+1, seq.AAMap[aa2]))
				break
			}
		}
	}
}
