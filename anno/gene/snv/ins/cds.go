package ins

import (
	"OpenAnno/db/transcript"
	"OpenAnno/seq"
	snv2 "OpenAnno/variant"
	"fmt"
)

func (a *GeneAnnoItem) AnnoInCDS(ins snv2.Snv, trans transcript.Transcript, regionIndex int, exonLen int) {
	region := trans.Regions[regionIndex]
	var pos int
	a.SetExon(region.ExonOrder)
	a.SetRegion("exonic")
	if trans.Strand == '+' {
		pos = exonLen + ins.Start - region.Start + 1
	} else {
		pos = exonLen + region.End - ins.Start
	}
	if trans.IsCmpl() {
		cdna, protein := trans.Cdna, trans.Protein
		newCdna := cdna.ChangeWithIns(pos, ins.Alt)
		newProtein := newCdna.Translate(ins.Chrom == "MT")

		for i := 0; i < cdna.Len(); i++ {
			if cdna.Base(i) != newCdna.Base(i) {
				a.SetAAChange(fmt.Sprintf("c.%dins%s", i, newCdna.SubSeq(i, ins.Alt.Len())))
				break
			}
		}
		lenL, lenR := 0, 0
		for i := 0; i < protein.Len(); i++ {
			if protein.Base(i) != newProtein.Base(i) {
				break
			}
			lenL++
		}
		if lenL < protein.Len() {
			if ins.Alt.Len()%3 == 0 {
				if newProtein.Find('*') < 0 {
					a.SetEvent("ins_nonframeshift_stoploss")
				} else if newProtein.Find('*') < protein.Len()-1 {
					a.SetEvent("ins_nonframeshift_stopgain")
				} else if protein.Base(0) != newProtein.Base(0) && protein.Base(0) == 'M' {
					a.SetEvent("ins_nonframeshift_startloss")
				} else {
					a.SetEvent("ins_nonframeshift")
				}
				for i := protein.Len() - 1; i >= 0; i-- {
					if protein.Base(i) != newProtein.Base(i) {
						break
					}
					lenR++
				}
				if lenL+lenR > protein.Len() {
					lenR = protein.Len() - lenL
				}
				start, end, varEnd := lenL, protein.Len()-lenR, newProtein.Len()-lenR
				altAa := newProtein.SubSeq(start, varEnd-start)
				if start == end {
					a.SetAAChange(fmt.Sprintf(
						"p.%s%d_%s%dins%s",
						seq.AAMap[protein.Base(start-1)],
						start,
						seq.AAMap[protein.Base(start)],
						start+1,
						altAa.ProteinOne2Tree(),
					))
				} else {
					refAa := protein.SubSeq(start, end-start)
					a.SetAAChange(fmt.Sprintf("p.%s%ddelins%s", refAa.ProteinOne2Tree(), start, altAa.ProteinOne2Tree()))
				}
			} else {
				if newProtein.Find('*') < 0 {
					a.SetEvent("ins_frameshift_stoploss")
				} else if newProtein.Find('*') < protein.Len()-1 {
					a.SetEvent("ins_frameshift_stopgain")
				} else if protein.Base(0) != newProtein.Base(0) {
					a.SetEvent("ins_frameshift_startloss")
				} else {
					a.SetEvent("ins_frameshift")
				}
				start := lenL
				a.SetAAChange(fmt.Sprintf(
					"p.%s%d%sfs",
					seq.AAMap[protein.Base(start)],
					start+1,
					seq.AAMap[newProtein.Base(start)],
				))
			}
		}
	}
}
