package gene

import (
	"fmt"
	"open-anno/pkg"
	"strings"
)

func getCDSPosOfNAchange(trans pkg.Transcript, cdsLen, pos, cLen int, region pkg.Region) (string, int) {
	if pos >= trans.CdsStart && pos <= trans.CdsEnd {
		if trans.Strand == "+" {
			if region.Type == pkg.RType_CDS {
				return fmt.Sprintf("%d", cLen+(pos-region.Start+1)), -1
			} else {
				dist1, dist2 := pos-region.Start+1, region.End-pos+1
				if dist1 < dist2 {
					return fmt.Sprintf("%d+%d", cLen, dist1), dist1
				} else {
					return fmt.Sprintf("%d-%d", cLen+1, dist2), dist2
				}
			}
		} else {
			if region.Type == pkg.RType_CDS {
				return fmt.Sprintf("%d", cdsLen-(cLen+pos-region.Start)), -1
			} else {
				dist1, dist2 := pos-region.Start+1, region.End-pos+1
				if dist1 < dist2 {
					return fmt.Sprintf("%d-%d", cdsLen-cLen+1, dist1), dist1
				} else {
					return fmt.Sprintf("%d+%d", cdsLen-cLen, dist2), dist2
				}
			}
		}
	}
	return "", -1
}

func setDelAAChange(transAnno TransAnno, trans pkg.Transcript, cstart, cend int) TransAnno {
	cdna := trans.CDNA()
	ncdna := pkg.Delete(cdna, cstart, cend)
	if trans.Strand == "-" {
		cdna = pkg.RevComp(cdna)
		ncdna = pkg.RevComp(ncdna)
	}
	start := pkg.DifferenceSimple(cdna, ncdna)
	alt := cdna[start-1 : start+cend-cstart]
	if transAnno.NAChange == "" {
		if cstart == cend {
			transAnno.NAChange = fmt.Sprintf("c.%ddel%s", start, alt)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.%d_%ddel%s", start, start+cend-cstart, alt)
		}
	}
	protein := pkg.Translate(cdna, trans.Chrom == "MT" || trans.Chrom == "chrM")
	nprotein := pkg.Translate(ncdna, trans.Chrom == "MT" || trans.Chrom == "chrM")
	start, end1, end2 := pkg.Difference(protein, nprotein)
	aa1 := protein[start-1 : end1]
	aa2 := nprotein[start-1 : end2]
	if (len(cdna)-len(ncdna))%3 == 0 {
		transAnno.Event = "del_inframe"
		if len(aa2) == 0 {
			if len(aa1) == 1 {
				transAnno.AAChange = fmt.Sprintf("p.%s%ddel", pkg.AAName(aa1, AA_SHORT), start)
			} else {
				transAnno.AAChange = fmt.Sprintf(
					"p.%s%d_%s%ddel",
					pkg.AAName(aa1[0], AA_SHORT),
					start,
					pkg.AAName(aa1[len(aa1)-1], AA_SHORT),
					end1,
				)
			}
			if aa1[len(aa1)-1] == '*' {
				transAnno.AAChange += "ext*?"
				transAnno.Event += "_stoploss"
			}

		} else {
			if len(aa1) == 1 {
				transAnno.AAChange = fmt.Sprintf("p.%s%ddelins%s", pkg.AAName(aa1, AA_SHORT), start, pkg.AAName(aa2, AA_SHORT))
			} else {
				transAnno.AAChange = fmt.Sprintf(
					"p.%s%d_%s%ddelins%s",
					pkg.AAName(aa1[0], AA_SHORT),
					start,
					pkg.AAName(aa1[len(aa1)-1], AA_SHORT),
					end1,
					pkg.AAName(aa2, AA_SHORT))
			}
		}
	} else {
		if start < len(protein) {
			transAnno.Event = "del_frameshift"
			if len(aa2) == 0 {
				transAnno.AAChange = fmt.Sprintf("p.%s%dfs", pkg.AAName(aa1[0], AA_SHORT), start)
			} else {
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
	}
	if protein[0] == 'M' && len(nprotein) > 0 && nprotein[0] != 'M' {
		transAnno.Event += "_startloss"
	}
	if strings.Contains(transAnno.Region, "splic") {
		transAnno.Event += "_splicing"
	}
	if protein[len(protein)-1] == '*' && strings.IndexByte(nprotein, '*') == -1 {
		transAnno.Event += "_stoploss"
	}
	return transAnno
}

func AnnoDel(snv pkg.AnnoVariant, trans pkg.Transcript) TransAnno {
	// cStart, cEnd, region1, region2, isExonSplicing := getDelCLen(trans, snv)
	utrLen1, utrLen2 := trans.ULen()
	cdsLen := trans.CLen()
	region1, cLen1, uLen1 := trans.Region(snv.Start)
	region2, cLen2, uLen2 := trans.Region(snv.End)
	utrPosOfNAchange1 := getUTRPosOfNAchange(trans, utrLen1, utrLen2, snv.Start, uLen1, region1)
	utrPosOfNAchange2 := getUTRPosOfNAchange(trans, utrLen1, utrLen2, snv.End, uLen2, region2)
	cdsPosOfNAchange1, dist1 := getCDSPosOfNAchange(trans, cdsLen, snv.Start, cLen1, region1)
	cdsPosOfNAchange2, dist2 := getCDSPosOfNAchange(trans, cdsLen, snv.End, cLen2, region2)
	transAnno := NewTransAnno(trans, region1, region2)
	if !region1.Equal(region2) {
		if region2.Start-region1.End > 1 {
			transAnno.Region = "deletion"
		} else {
			if (region1.Type == pkg.RType_CDS && region2.Type == pkg.RType_INTRON) || (region2.Type == pkg.RType_CDS && region1.Type == pkg.RType_INTRON) {
				transAnno.Region = "exonic_splicing"
			}
		}
	}
	if (snv.Start < trans.CdsStart || snv.Start > trans.CdsEnd) && (snv.End < trans.CdsStart || snv.End > trans.CdsEnd) {
		// snv包含了整个编码区，即整个编码区被删除 SNV的起始终止位置都不在CDS区或CDS区之间的Intron中
		// 情况2: snv包含了整个编码区，即整个编码区被删除
		// ...+++,,,+++...
		//   |---------|
		//|---------------|
		// 情况2: snv发生在整个编码区左边
		// |-|...+++,,,+++...
		//   |-|
		// 情况3: snv发生在整个编码区左边
		// ...+++,,,+++...|-|
		//              |-|
		if utrPosOfNAchange1 == utrPosOfNAchange2 {
			transAnno.NAChange = fmt.Sprintf("c.%sdel", utrPosOfNAchange1)
		} else {
			if trans.Strand == "+" {
				transAnno.NAChange = fmt.Sprintf("c.%s_%sdel", utrPosOfNAchange1, utrPosOfNAchange2)
			} else {
				transAnno.NAChange = fmt.Sprintf("c.%s_%sdel", utrPosOfNAchange2, utrPosOfNAchange1)
			}
		}
		if snv.Start < trans.CdsStart && snv.End > trans.CdsEnd {
			// 情况1时
			transAnno.Region = "transcript"
			transAnno.Region2 = "transcript"
			transAnno.Event = "deletion"
		}
	} else if trans.CdsStart > snv.Start && trans.CdsStart <= snv.End && trans.CdsEnd >= snv.End {
		// ...+++,,,+++...
		//  |--|
		//|----|
		if trans.Strand == "+" {
			transAnno.NAChange = fmt.Sprintf("c.%s_%sdel", utrPosOfNAchange1, cdsPosOfNAchange2)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.%s_%sdel", cdsPosOfNAchange2, utrPosOfNAchange1)
		}
	} else if trans.CdsEnd >= snv.Start && trans.CdsEnd < snv.End && trans.CdsStart <= snv.Start {
		// ...,,,+++...
		//        |--|
		//        |----|
		if trans.Strand == "+" {
			transAnno.NAChange = fmt.Sprintf("c.%s_%sdel", cdsPosOfNAchange1, utrPosOfNAchange2)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.%s_%sdel", utrPosOfNAchange2, cdsPosOfNAchange1)
		}

	} else {
		// 变异跨region
		if !region1.Equal(region2) {
			if cdsPosOfNAchange1 == cdsPosOfNAchange2 {
				transAnno.NAChange = fmt.Sprintf("c.%sdel", cdsPosOfNAchange1)
			} else {
				if trans.Strand == "+" {
					transAnno.NAChange = fmt.Sprintf("c.%s_%sdel", cdsPosOfNAchange1, cdsPosOfNAchange2)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%s_%sdel", cdsPosOfNAchange2, cdsPosOfNAchange1)
				}
			}
		} else {
			if region1.Type == pkg.RType_INTRON {
				if cdsPosOfNAchange1 == cdsPosOfNAchange2 {
					transAnno.NAChange = fmt.Sprintf("c.%sdel", cdsPosOfNAchange1)
				} else {
					if trans.Strand == "+" {
						transAnno.NAChange = fmt.Sprintf("c.%s_%sdel", cdsPosOfNAchange1, cdsPosOfNAchange2)
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%s_%sdel", cdsPosOfNAchange2, cdsPosOfNAchange1)
					}
				}
				if dist1 > 0 && dist2 > 0 && pkg.Min(dist1, dist2) <= 2 {
					transAnno.Event = "splicing"
					transAnno.Region = "splicing"
				}
			} else {
				cstart, cend := cLen1, cLen2
				cstart += snv.Start - region1.Start + 1
				cend += snv.End - region2.Start + 1
				transAnno = setDelAAChange(transAnno, trans, cstart, cend)
			}
		}
	}
	return transAnno
}

func AnnoUnkDel(snv pkg.AnnoVariant, trans pkg.Transcript) TransAnno {
	transAnno := NewTransAnno(trans)
	transAnno.Region2 = "ncRNA"
	if snv.Start < trans.TxStart {
		if snv.End > trans.TxEnd {
			if trans.Strand == "+" {
				transAnno.NAChange = fmt.Sprintf("n.-%d_+%ddel", trans.TxStart-snv.Start, snv.End-trans.TxEnd)
			} else {
				transAnno.NAChange = fmt.Sprintf("n.-%d_+%ddel", snv.End-trans.TxEnd, trans.TxStart-snv.Start)
			}
		} else {
			if trans.Strand == "+" {
				transAnno.NAChange = fmt.Sprintf("n.-%d_%ddel", trans.TxStart-snv.Start, snv.End-trans.TxStart+1)
			} else {
				transAnno.NAChange = fmt.Sprintf("n.%d_+%ddel", trans.TxEnd-trans.TxEnd+1, trans.TxStart-snv.Start)
			}
		}
	} else {
		if snv.End > trans.TxEnd {
			if trans.Strand == "+" {
				transAnno.NAChange = fmt.Sprintf("n.%d_+%ddel", trans.TxStart-snv.Start+1, snv.End-trans.TxEnd)
			} else {
				transAnno.NAChange = fmt.Sprintf("n.-%d_%ddel", snv.End-trans.TxEnd, trans.TxEnd-snv.Start+1)
			}
		} else {
			nstart := snv.Start - trans.TxStart + 1
			nend := snv.End - trans.TxStart + 1
			dna := trans.DNA()
			ndna := pkg.Delete(dna, nstart, nend)
			if trans.Strand == "-" {
				dna = pkg.RevComp(dna)
				ndna = pkg.RevComp(ndna)
			}
			start := pkg.DifferenceSimple(dna, ndna)
			// alt := dna[start-1 : start+nend-nstart]
			if nstart == nend {
				// transAnno.NAChange = fmt.Sprintf("n.%ddel%s", start, alt)
				transAnno.NAChange = fmt.Sprintf("n.%ddel", start)
			} else {
				// transAnno.NAChange = fmt.Sprintf("n.%d_%ddel%s", start, start+nend-nstart, alt)
				transAnno.NAChange = fmt.Sprintf("n.%d_%ddel", start, start+nend-nstart)
			}
		}
	}
	return transAnno
}
