package gene

import (
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"open-anno/pkg/io/refgene"
	"open-anno/pkg/seq"
	"strings"
)

func findInsRegion(regions refgene.Regions, strand string, snv io.Variant) (refgene.Region, int) {
	var cLen int
	for _, region := range regions {
		if (strand == "+" && snv.Start >= region.Start-1 && snv.End < region.End) ||
			(strand == "-" && snv.Start >= region.Start && snv.End <= region.End) {
			return region, cLen
		}
		if region.Type == refgene.RType_CDS {
			cLen += region.End - region.Start + 1
		}
	}
	return refgene.Region{}, cLen
}

func AnnoIns(snv io.Variant, trans refgene.Transcript, aashort bool) SnvGeneBased {
	region, cLen := findInsRegion(trans.Regions, trans.Strand, snv)
	anno := NewSnvGeneBased(trans, region)
	if region.End < trans.CdsStart {
		if trans.Strand == "+" {
			anno.NAChange = fmt.Sprintf("c.-%d_-%dins%s", trans.CdsStart-snv.End, trans.CdsStart-snv.End-1, snv.Alt)
		} else {
			anno.NAChange = fmt.Sprintf("c.+%d_+%dins%s", trans.CdsStart-snv.End-1, trans.CdsStart-snv.End, seq.RevComp(snv.Alt))
		}
	} else if region.Start > trans.CdsEnd {
		if trans.Strand == "+" {
			anno.NAChange = fmt.Sprintf("c.+%d_+%dins%s", snv.Start-trans.CdsEnd, snv.Start-trans.CdsEnd+1, snv.Alt)
		} else {
			anno.NAChange = fmt.Sprintf("c.-%d_-%dins%s", snv.Start-trans.CdsEnd+1, snv.Start-trans.CdsEnd+1, seq.RevComp(snv.Alt))
		}
	} else {
		if region.Type == refgene.RType_INTRON {
			dist1s, dist2s := snv.Start-region.Start+1, region.End-snv.End+1
			dist1e, dist2e := dist1s+1, dist2s-1
			dist1, dist2 := dist1e, dist2s
			if pkg.Min(dist1, dist2) <= 2 {
				anno.Event = "splicing"
			}
			if trans.Strand == "+" {
				if dist1 <= dist2 {
					if dist1s == 0 {
						anno.NAChange = fmt.Sprintf("c.%d_%d+%dins%s", cLen, cLen, dist1e, snv.Alt)
					} else {
						anno.NAChange = fmt.Sprintf("c.%d+%d_%d+%dins%s", cLen, dist1s, cLen, dist1e, snv.Alt)
					}
				} else {
					if dist2e == 0 {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%dins%s", cLen+1, cLen+1, dist2s, snv.Alt)
					} else {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%d-%dins%s", cLen+1, dist2s, cLen+1, dist2e, snv.Alt)
					}
				}
			} else {
				if dist1 <= dist2 {
					nclen := trans.CLen() - cLen + 1
					if dist1s == 0 {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%dins%s", nclen, dist1e, nclen, seq.RevComp(snv.Alt))
					} else {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%d-%dins%s", nclen, dist1e, nclen, dist1s, seq.RevComp(snv.Alt))
					}
				} else {
					nclen := trans.CLen() - cLen
					if dist2e == 0 {
						anno.NAChange = fmt.Sprintf("c.%d_%d+%dins%s", nclen, nclen, dist2s, seq.RevComp(snv.Alt))
					} else {
						anno.NAChange = fmt.Sprintf("c.%d+%d_%d+%dins%s", nclen, dist2e, nclen, dist2s, seq.RevComp(snv.Alt))
					}
				}
			}
		} else {
			pos := cLen + snv.Start - region.Start + 1
			cdna := trans.CDNA()
			ncdna := seq.Insert(cdna, pos, snv.Alt)
			if trans.Strand == "-" {
				cdna = seq.RevComp(cdna)
				ncdna = seq.RevComp(ncdna)
			}
			protein := seq.Translate(cdna, trans.Chrom == "MT")
			nprotein := seq.Translate(ncdna, trans.Chrom == "MT")
			start := seq.DifferenceSimple(cdna, ncdna)
			alt := ncdna[start-1 : start+len(snv.Alt)-1]
			unit := seq.DupUnit(alt)
			pre := cdna[start-len(unit)-1 : start-1]
			if unit == pre {
				anno.NAChange = fmt.Sprintf("c.%ddup%s", start-1, unit)
			} else {
				anno.NAChange = fmt.Sprintf("c.%d_%dins%s", start-1, start, alt)
			}
			start, end1, end2 := seq.Difference(protein, nprotein)
			aa1 := protein[start-1 : end1]
			aa2 := nprotein[start-1 : end2]
			if len(snv.Alt)%3 == 0 {
				anno.Event = "ins_nonframeshift"
				if len(aa1) == 0 {
					anno.AAChange = fmt.Sprintf(
						"p.%s%d_%s%dins%s",
						seq.AAName(protein[start-2], aashort),
						start-1,
						seq.AAName(protein[start-1], aashort),
						start,
						seq.AAName(aa2, aashort),
					)
				} else if len(aa1) == 1 {
					anno.AAChange = fmt.Sprintf("p.%s%ddelins%s", seq.AAName(aa1, aashort), start-1, seq.AAName(aa2, aashort))
				} else {
					anno.AAChange = fmt.Sprintf(
						"p.%s%d_%s%ddelins%s",
						seq.AAName(aa1[0], aashort),
						start-1,
						seq.AAName(aa1[len(aa1)-1], aashort),
						end1-1,
						seq.AAName(aa2, aashort))
				}
			} else {
				if start < len(protein) {
					anno.Event = "ins_frameshift"
					var fs string
					fsi := strings.IndexByte(nprotein[start-1:], '*')
					if fsi == -1 {
						fs = "?"
					}
					if fsi != 0 {
						fs = fmt.Sprintf("%d", fsi+1)
					}
					anno.AAChange = fmt.Sprintf("p.%s%d%sfs*%s", seq.AAName(aa1[0], aashort), start, seq.AAName(aa2[0], aashort), fs)
				}
			}
		}
	}
	return anno
}
