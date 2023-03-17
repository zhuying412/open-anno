package gene

import (
	"fmt"
	"open-anno/pkg"
	"strings"
)

func findInsRegion(regions pkg.Regions, strand string, snv pkg.AnnoVariant) (pkg.Region, int) {
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

func AnnoIns(snv pkg.AnnoVariant, trans pkg.Transcript) TransAnno {
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
			if start-1 == len(cdna) {

			}
			alt := ncdna[start-1 : start+len(snv.Alt)-1]
			if start == 1 {
				// 在CDNA第一个碱基前插入，不影响起始密码子
				transAnno.NAChange = fmt.Sprintf("c.-1_1ins%s", alt)
				if trans.HasUTR5() {
					transAnno.Region = "UTR5"
					transAnno.Region2 = "UTR5"
				} else {
					transAnno.Region = "UpStream"
					transAnno.Region2 = "UpStream"
				}

			} else if start-1 == len(cdna) {
				// 在CDNA最后碱基后插入，不影响整条蛋白序列
				transAnno.NAChange = fmt.Sprintf("c.%d_%d+1ins%s", len(cdna), len(cdna), alt)
				if trans.HasUTR3() {
					transAnno.Region = "UTR3"
					transAnno.Region2 = "UTR3"
				} else {
					transAnno.Region = "DownStream"
					transAnno.Region2 = "DownStream"
				}
			} else {
				unit := pkg.DupUnit(alt)
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
						unit = pkg.DupUnit(aa2)
						if start > len(unit) && unit == cdna[start-len(unit)-1:start-1] {
							// 比较重复单元和其之前的碱基
							if len(unit) == 1 {
								transAnno.NAChange = fmt.Sprintf("c.%s%ddup", unit, start-1)
							} else {
								transAnno.NAChange = fmt.Sprintf(
									"c.%s%d_%s%ddup",
									pkg.AAName(protein[start-len(unit)-1], AA_SHORT),
									start-len(unit),
									pkg.AAName(protein[start-2], AA_SHORT),
									start-1,
								)
							}
						} else {
							transAnno.AAChange = fmt.Sprintf(
								"p.%s%d_%s%dins%s",
								pkg.AAName(protein[start-2], AA_SHORT),
								start-1,
								pkg.AAName(protein[start-1], AA_SHORT),
								start,
								pkg.AAName(aa2, AA_SHORT),
							)
						}
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
						if aa2[0] == '*' {
							transAnno.AAChange = fmt.Sprintf("p.%s%d*", pkg.AAName(aa1[0], AA_SHORT), start)
						} else {
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
				if protein[0] == 'M' && nprotein[0] != 'M' {
					transAnno.Event += "_startloss"
				}
				if protein[len(protein)-1] == '*' && strings.IndexByte(nprotein, '*') == -1 {
					transAnno.Event += "_stoploss"
				}
			}
		}
	}
	return transAnno
}
