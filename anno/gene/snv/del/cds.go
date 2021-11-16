package del

import (
	"OpenAnno/db/transcript"
	"OpenAnno/seq"
	"OpenAnno/variant"
	"fmt"
)

func GetExonLen(lenL int, lenR int, start int, end int, regions transcript.Regions, strand byte) (int, int) {
	for _, region := range regions {
		if region.IsCDS() {
			if region.Start <= start && start <= region.End {
				if strand == '+' {
					lenL += start - region.Start
				} else {
					lenR += start - region.Start
				}
			}
			if region.Start <= end && end <= region.End {
				if strand == '+' {
					lenR += region.End - end
				} else {
					lenL += region.End - end
				}
			}
		}
	}
	return lenL, lenR
}

func (a *GeneAnnoItem) AnnoInCDS(del variant.Snv, trans transcript.Transcript, lenL int, lenR int, regions transcript.Regions) {
	lenL, lenR = GetExonLen(lenL, lenR, del.Start, del.End, regions, trans.Strand)
	cdna, protein := trans.Cdna, trans.Protein
	newCdna := cdna.ChangeWithDel(lenL, lenR)
	newProtein := newCdna.Translate(del.Chrom == "MT")
	lenl, lenr := 0, 0
	for i := 0; i < lenL+lenR; i++ {
		if cdna[i] != newCdna[i] {
			break
		}
		lenl++
	}
	for i := lenL + lenR - 1; i >= 0; i-- {
		if cdna[i] != newCdna[i] {
			break
		}
		lenr++
	}
	if lenl+lenr > lenL+lenR {
		lenr = lenL + lenR - lenl
	}
	a.SetAAChange(fmt.Sprintf("c.%d_%ddel", lenl+1, lenL+lenR-lenr))
	lenDel := lenL + lenR - lenl - lenr
	lenl, lenr, lenp, lenvp := 0, 0, protein.Len(), newProtein.Len()
	for i := 0; i < lenvp; i++ {
		if protein[i] == newProtein[i] {
			break
		}
		lenl++
	}
	if lenDel%3 == 0 {
		for i := lenvp - 1; i >= 0; i-- {
			if protein[i] != newProtein[i] {
				break
			}
			lenr++
		}
		if lenl+lenr > lenvp {
			lenr = lenvp - lenl
		}
		start, end1, end2 := lenl+1, lenp-lenr, lenvp-lenr
		if newProtein.Find('*') < 0 {
			a.SetEvent("del_nonframeshift_stoploss")
		} else if newProtein.Find('*') < lenvp-1 {
			a.SetEvent("del_nonframeshift_stopgain")
		} else if protein.Base(0) != protein.Base(0) {
			a.SetEvent("del_nonframeshift_startloss")
		} else {
			a.SetEvent("del_nonframeshift")
		}
		if start == end2+1 {
			if start == end1 {
				a.SetAAChange(fmt.Sprintf("p.%s%ddel", seq.AAMap[protein.Base(start-1)], start))
			} else {
				a.SetAAChange(fmt.Sprintf("p.%s%d_%s%ddel",
					seq.AAMap[protein.Base(start-1)],
					start,
					seq.AAMap[protein.Base(end1-1)],
					end1))
			}
		} else {
			a.SetAAChange(fmt.Sprintf("p.%s%d_%s%dinsdel%s",
				seq.AAMap[protein.Base(start-1)],
				start,
				seq.AAMap[protein.Base(end1-1)],
				end1,
				newProtein.SubSeq(start-1, end2-start+1).ProteinOne2Tree()))
		}
	} else {
		start := lenl + 1
		if start > newProtein.Len() {
			a.SetEvent("del_nonframeshift_stoploss")
			a.SetAAChange(fmt.Sprintf("p.%s%d_%s%s%ddel", seq.AAMap[protein.Base(start-1)], start, seq.AAMap[protein[lenp-1]], lenp))
		} else {
			if newProtein.Find('*') < 0 {
				a.SetEvent("del_frameshift_stoploss")
			} else if newProtein.Find('*') < lenvp-1 {
				a.SetEvent("del_frameshift_stopgain")
			} else if protein.Base(0) != protein.Base(0) {
				a.SetEvent("del_nonframeshift_startloss")
			} else {
				a.SetEvent("del_frameshift")
			}
			a.SetAAChange(fmt.Sprintf("p.%s%d%sfs",
				seq.AAMap[protein.Base(start-1)],
				start,
				seq.AAMap[newProtein.Base(start-1)]))
		}
	}
	a.SetRegion("exonic")
	for _, region := range regions {
		if region.IsCDS() {
			a.SetExon(region.ExonOrder)
			break
		}
	}
}
