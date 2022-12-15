package gene

import (
	"fmt"
	"open-anno/anno"
	"open-anno/pkg"
	"strings"
)

func findInsRegion(regions pkg.Regions, strand string, snv anno.Variant) (pkg.Region, int) {
	var cLen int
	for _, region := range regions {

		if (strand == "+" && snv.Start >= region.Start-1 && snv.End < region.End) ||
			(strand == "-" && snv.Start >= region.Start && snv.End <= region.End) {
			return region, cLen
		}
		if region.Type == pkg.RType_CDS {
			cLen += region.End - region.Start + 1
		}
	}
	return pkg.Region{}, cLen
}

func AnnoIns(snv anno.Variant, trans pkg.Transcript) TransAnno {
	region, cLen := findInsRegion(trans.Regions, trans.Strand, snv)
	transAnno := NewTransAnno(trans, region)
	if region.End < trans.CdsStart {
		if trans.Strand == "+" {
			if trans.CdsStart-snv.End > 1 {
				transAnno.NAChange = fmt.Sprintf("c.-%d_-%dins%s", trans.CdsStart-snv.End, trans.CdsStart-snv.End-1, snv.Alt)
			} else {
				transAnno.NAChange = fmt.Sprintf("c.-%d_1ins%s", trans.CdsStart-snv.End, snv.Alt)
			}
		} else {
			if trans.CdsStart-snv.End > 1 {
				transAnno.NAChange = fmt.Sprintf("c.+%d_+%dins%s", trans.CdsStart-snv.End-1, trans.CdsStart-snv.End, pkg.RevComp(snv.Alt))
			} else {
				transAnno.NAChange = fmt.Sprintf("c.%d_+%dins%s", trans.CLen(), trans.CdsStart-snv.End, pkg.RevComp(snv.Alt))
			}
		}
	} else if region.Start > trans.CdsEnd {
		if trans.Strand == "+" {
			if snv.Start-trans.CdsEnd > 0 {
				transAnno.NAChange = fmt.Sprintf("c.+%d_+%dins%s", snv.Start-trans.CdsEnd, snv.Start-trans.CdsEnd+1, snv.Alt)
			} else {
				transAnno.NAChange = fmt.Sprintf("c.%d_+%dins%s", trans.CLen(), snv.Start-trans.CdsEnd+1, snv.Alt)
			}
		} else {
			if snv.Start-trans.CdsEnd > 0 {
				transAnno.NAChange = fmt.Sprintf("c.-%d_-%dins%s", snv.Start-trans.CdsEnd+1, snv.Start-trans.CdsEnd, pkg.RevComp(snv.Alt))
			} else {
				transAnno.NAChange = fmt.Sprintf("c.-%d_1ins%s", snv.Start-trans.CdsEnd+1, pkg.RevComp(snv.Alt))
			}
		}
	} else {
		if region.Type == pkg.RType_INTRON {
			dist1s, dist2s := snv.Start-region.Start+1, region.End-snv.End+1
			dist1e, dist2e := dist1s+1, dist2s-1
			dist1, dist2 := dist1e, dist2s
			if pkg.Min(dist1, dist2) <= 2 {
				transAnno.Event = "splicing"
				transAnno.Region = "splicing"
			}
			if trans.Strand == "+" {
				if dist1 <= dist2 {
					if dist1s == 0 {
						transAnno.NAChange = fmt.Sprintf("c.%d_%d+%dins%s", cLen, cLen, dist1e, snv.Alt)
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%d+%d_%d+%dins%s", cLen, dist1s, cLen, dist1e, snv.Alt)
					}
				} else {
					if dist2e == 0 {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%dins%s", cLen+1, cLen+1, dist2s, snv.Alt)
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%d-%dins%s", cLen+1, dist2s, cLen+1, dist2e, snv.Alt)
					}
				}
			} else {
				if dist1 <= dist2 {
					nclen := trans.CLen() - cLen + 1
					if dist1s == 0 {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%dins%s", nclen, dist1e, nclen, pkg.RevComp(snv.Alt))
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%d-%dins%s", nclen, dist1e, nclen, dist1s, pkg.RevComp(snv.Alt))
					}
				} else {
					nclen := trans.CLen() - cLen
					if dist2e == 0 {
						transAnno.NAChange = fmt.Sprintf("c.%d_%d+%dins%s", nclen, nclen, dist2s, pkg.RevComp(snv.Alt))
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%d+%d_%d+%dins%s", nclen, dist2e, nclen, dist2s, pkg.RevComp(snv.Alt))
					}
				}
			}
		} else {
			pos := cLen + snv.Start - region.Start + 1
			cdna := trans.CDNA()
			ncdna := pkg.Insert(cdna, pos, snv.Alt)
			if trans.Strand == "-" {
				cdna = pkg.RevComp(cdna)
				ncdna = pkg.RevComp(ncdna)
			}
			protein := pkg.Translate(cdna, trans.Chrom == "MT")
			nprotein := pkg.Translate(ncdna, trans.Chrom == "MT")
			start := pkg.DifferenceSimple(cdna, ncdna)
			alt := ncdna[start-1 : start+len(snv.Alt)-1]
			unit := pkg.DupUnit(alt)
			if start == 1 {
				// 在CDNA第一个碱基前插入，不影响起始密码子
				transAnno.NAChange = fmt.Sprintf("c.-1_1ins%s", alt)
			} else {
				if start > len(unit) && unit == cdna[start-len(unit)-1:start-1] {
					// 比较重复单元和其之前的碱基
					transAnno.NAChange = fmt.Sprintf("c.%ddup%s", start-1, unit)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d_%dins%s", start-1, start, alt)
				}
				start, end1, end2 := pkg.Difference(protein, nprotein)
				aa1 := protein[start-1 : end1]
				aa2 := nprotein[start-1 : end2]
				if len(snv.Alt)%3 == 0 {
					transAnno.Event = "ins_inframe"
					if len(aa1) == 0 {
						transAnno.AAChange = fmt.Sprintf(
							"p.%s%d_%s%dins%s",
							pkg.AAName(protein[start-2], AA_SHORT),
							start-1,
							pkg.AAName(protein[start-1], AA_SHORT),
							start,
							pkg.AAName(aa2, AA_SHORT),
						)
					} else if len(aa1) == 1 {
						transAnno.AAChange = fmt.Sprintf("p.%s%ddelins%s", pkg.AAName(aa1, AA_SHORT), start-1, pkg.AAName(aa2, AA_SHORT))
					} else {
						transAnno.AAChange = fmt.Sprintf(
							"p.%s%d_%s%ddelins%s",
							pkg.AAName(aa1[0], AA_SHORT),
							start-1,
							pkg.AAName(aa1[len(aa1)-1], AA_SHORT),
							end1-1,
							pkg.AAName(aa2, AA_SHORT))
					}
				} else {
					if start < len(protein) {
						transAnno.Event = "ins_frameshift"
						var fs string
						fsi := strings.IndexByte(nprotein[start-1:], '*')
						if fsi == -1 {
							fs = "?"
						}
						if fsi != 0 {
							fs = fmt.Sprintf("%d", fsi+1)
						}
						transAnno.AAChange = fmt.Sprintf("p.%s%d%sfs*%s", pkg.AAName(aa1[0], AA_SHORT), start, pkg.AAName(aa2[0], AA_SHORT), fs)
					}
				}
			}
		}
	}
	return transAnno
}
