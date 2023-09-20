package gene

import (
	"fmt"
	"open-anno/pkg"
	"strings"
)

func setInsAAChange(transAnno TransAnno, trans pkg.Transcript, snv pkg.AnnoVariant, cPos int) TransAnno {
	// pos := cLen + snv.Start - region.Start + 1
	cdna := trans.CDNA()
	ncdna := pkg.Insert(cdna, cPos, snv.Alt)
	if trans.Strand == "-" {
		cdna = pkg.RevComp(cdna)
		ncdna = pkg.RevComp(ncdna)
	}
	protein := pkg.Translate(cdna, trans.Chrom == "MT" || trans.Chrom == "chrM")
	nprotein := pkg.Translate(ncdna, trans.Chrom == "MT" || trans.Chrom == "chrM")
	start := pkg.DifferenceSimple(cdna, ncdna)
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
			if len(unit) == 1 {
				transAnno.NAChange = fmt.Sprintf("c.%ddup%s", start-1, unit)
			} else {
				transAnno.NAChange = fmt.Sprintf("c.%d_%ddup%s", start-len(unit), start-1, unit)
			}

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
				if start > len(unit) && unit == protein[start-len(unit)-1:start-1] {
					// 比较重复单元和其之前的碱基
					if len(unit) == 1 {
						transAnno.AAChange = fmt.Sprintf("p.%s%ddup", unit, start-1)
					} else {
						transAnno.AAChange = fmt.Sprintf(
							"p.%s%d_%s%ddup",
							pkg.AAName(protein[start-len(unit)-1], AA_SHORT),
							start-len(unit),
							pkg.AAName(protein[start-2], AA_SHORT),
							start-1,
						)
					}
				} else {
					if start == 1 {
						transAnno.AAChange = fmt.Sprintf(
							"p.%s%d-1_%s%dins%s",
							pkg.AAName(protein[start-1], AA_SHORT),
							start,
							pkg.AAName(protein[start-1], AA_SHORT),
							start,
							pkg.AAName(aa2, AA_SHORT),
						)
					} else {
						if start > len(protein) {
							transAnno.AAChange = fmt.Sprintf(
								"p.%s%d_%s%d+1ins%s",
								pkg.AAName(protein[start-2], AA_SHORT),
								start-1,
								pkg.AAName(protein[start-2], AA_SHORT),
								start-1,
								pkg.AAName(aa2, AA_SHORT),
							)
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
					}
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
					} else {
						if fsi != 0 {
							fs = fmt.Sprintf("%d", fsi+1)
						}
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
	return transAnno
}

func AnnoIns(snv pkg.AnnoVariant, trans pkg.Transcript) TransAnno {
	utrLen1, utrLen2 := trans.ULen()
	cdsLen := trans.CLen()
	region1, cLen1, uLen1 := trans.Region(snv.Start)
	region2, cLen2, uLen2 := trans.Region(snv.Start + 1)
	utrPosOfNAchange1 := getUTRPosOfNAchange(trans, utrLen1, utrLen2, snv.Start, uLen1, region1)
	utrPosOfNAchange2 := getUTRPosOfNAchange(trans, utrLen1, utrLen2, snv.Start+1, uLen2, region2)
	cdsPosOfNAchange1, dist1 := getCDSPosOfNAchange(trans, cdsLen, snv.Start, cLen1, region1)
	cdsPosOfNAchange2, dist2 := getCDSPosOfNAchange(trans, cdsLen, snv.Start+1, cLen2, region2)
	var transAnno TransAnno
	if trans.Strand == "+" {
		transAnno = NewTransAnno(trans, region2)
	} else {
		transAnno = NewTransAnno(trans, region1)
	}
	if !region1.Equal(region2) {
		if !region1.Exists() || !region2.Exists() || region2.End < trans.CdsStart || region1.Start > trans.CdsEnd {
			if trans.Strand == "+" {
				transAnno.NAChange = fmt.Sprintf("c.%s_%sins%s", utrPosOfNAchange1, utrPosOfNAchange2, snv.Alt)
			} else {
				transAnno.NAChange = fmt.Sprintf("c.%s_%sins%s", utrPosOfNAchange2, utrPosOfNAchange1, pkg.RevComp(snv.Alt))
			}
		} else if region1.Start >= trans.CdsStart && region2.End <= trans.CdsEnd {
			if trans.Strand == "+" {
				if region2.Type == pkg.RType_CDS {
					transAnno = setInsAAChange(transAnno, trans, snv, cLen2)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%s_%sins%s", cdsPosOfNAchange1, cdsPosOfNAchange2, snv.Alt)
					transAnno.Event = "splicing"
					transAnno.Region = "splicing"
				}
			} else {
				if region1.Type == pkg.RType_CDS {
					transAnno = setInsAAChange(transAnno, trans, snv, cLen2)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%s_%sins%s", cdsPosOfNAchange2, cdsPosOfNAchange1, pkg.RevComp(snv.Alt))
					transAnno.Event = "splicing"
					transAnno.Region = "splicing"
				}
			}

		} else {
			if trans.Strand == "+" {
				if region1.End < trans.CdsStart {
					transAnno = setInsAAChange(transAnno, trans, snv, cLen2)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%s_%sins%s", cdsPosOfNAchange1, utrPosOfNAchange2, snv.Alt)
				}
			} else {
				if region1.End < trans.CdsEnd {
					transAnno = setInsAAChange(transAnno, trans, snv, cLen2)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%s_%sins%s", cdsPosOfNAchange2, utrPosOfNAchange1, pkg.RevComp(snv.Alt))
				}
			}
		}
	} else {
		if region1.End < trans.CdsStart || region1.Start > trans.CdsEnd {
			if trans.Strand == "+" {
				transAnno.NAChange = fmt.Sprintf("c.%s_%sins%s", utrPosOfNAchange1, utrPosOfNAchange2, snv.Alt)
			} else {
				transAnno.NAChange = fmt.Sprintf("c.%s_%sins%s", utrPosOfNAchange2, utrPosOfNAchange1, pkg.RevComp(snv.Alt))
			}
		} else {
			if region1.Type == pkg.RType_CDS {
				transAnno = setInsAAChange(transAnno, trans, snv, cLen1+snv.Start-region1.Start+1)
			} else {
				if trans.Strand == "+" {
					transAnno.NAChange = fmt.Sprintf("c.%s_%sins%s", cdsPosOfNAchange1, cdsPosOfNAchange2, snv.Alt)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%s_%sins%s", cdsPosOfNAchange2, cdsPosOfNAchange1, pkg.RevComp(snv.Alt))
				}
				if pkg.Min(dist1, dist2) <= 2 {
					transAnno.Event = "splicing"
					transAnno.Region = "splicing"
				}
			}
		}
	}
	return transAnno
}

func AnnoUnkIns(snv pkg.AnnoVariant, trans pkg.Transcript) TransAnno {
	transAnno := NewTransAnno(trans)
	transAnno.Region2 = "ncRNA"
	pos := snv.Start - trans.TxStart + 1
	dna := trans.DNA()
	ndna := pkg.Insert(dna, pos, snv.Alt)
	if trans.Strand == "-" {
		dna = pkg.RevComp(dna)
		ndna = pkg.RevComp(ndna)
	}
	start := pkg.DifferenceSimple(dna, ndna)
	alt := ndna[start-1 : start+len(snv.Alt)-1]
	if start == 1 {
		// 在CDNA第一个碱基前插入，不影响起始密码子
		transAnno.NAChange = fmt.Sprintf("n.-1_1ins%s", alt)
	} else if start-1 == len(dna) {
		// 在CDNA最后碱基后插入，不影响整条蛋白序列
		transAnno.NAChange = fmt.Sprintf("n.%d_%d+1ins%s", len(dna), len(dna), alt)
	} else {
		unit := pkg.DupUnit(alt)
		if start > len(unit) && unit == dna[start-len(unit)-1:start-1] {
			// 比较重复单元和其之前的碱基
			transAnno.NAChange = fmt.Sprintf("n.%ddup%s", start-1, unit)
		} else {
			transAnno.NAChange = fmt.Sprintf("n.%d_%dins%s", start-1, start, alt)
		}
	}
	return transAnno
}
