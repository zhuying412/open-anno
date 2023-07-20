package gene

import (
	"fmt"
	"open-anno/pkg"
	"strings"
)

func setSubAAChange(transAnno TransAnno, trans pkg.Transcript, cstart int, cend int, alt string) TransAnno {
	cdna := trans.CDNA()
	ncdna := pkg.Substitute2(cdna, cstart, cend, alt)
	if trans.Strand == "-" {
		alt = pkg.RevComp(alt)
		cdna = pkg.RevComp(cdna)
		ncdna = pkg.RevComp(ncdna)
	}
	start := pkg.DifferenceSimple(cdna, ncdna)
	if transAnno.NAChange == "" {
		if cstart == cend {
			transAnno.NAChange = fmt.Sprintf("c.%ddelins%s", start, alt)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", start, start+cend-cstart, alt)
		}

	}
	protein := pkg.Translate(cdna, trans.Chrom == "MT" || trans.Chrom == "chrM")
	nprotein := pkg.Translate(ncdna, trans.Chrom == "MT" || trans.Chrom == "chrM")
	start, end1, end2 := pkg.Difference(protein, nprotein)
	aa1 := protein[start-1 : end1]
	aa2 := nprotein[start-1 : end2]
	if (len(cdna)-len(ncdna))%3 == 0 {
		transAnno.Event = "sub_inframe"
		if len(aa2) == 0 {
			if len(aa1) == 1 {
				transAnno.AAChange = fmt.Sprintf("p.%s%ddel", pkg.AAName(aa1, AA_SHORT), start)
			} else if len(aa1) > 1 {
				transAnno.AAChange = fmt.Sprintf(
					"p.%s%d_%s%ddel",
					pkg.AAName(aa1[0], AA_SHORT),
					start,
					pkg.AAName(aa1[len(aa1)-1], AA_SHORT),
					end1,
				)
			}
		} else {
			if len(aa1) == 1 {
				transAnno.AAChange = fmt.Sprintf("p.%s%ddelins%s", pkg.AAName(aa1, AA_SHORT), start, pkg.AAName(aa2, AA_SHORT))
			} else if len(aa1) > 1 {
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
			transAnno.Event = "sub_frameshift"
			if aa2[0] == '*' {
				transAnno.AAChange = fmt.Sprintf("p.%s%dfs", pkg.AAName(aa1[0], AA_SHORT), start)
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
	if protein[0] != nprotein[0] && protein[0] == 'M' {
		transAnno.Event += "_startloss"
	}
	if strings.Contains(transAnno.Region, "splic") {
		transAnno.Event += "_splicing"
	}
	return transAnno
}

func AnnoSub(snv pkg.AnnoVariant, trans pkg.Transcript) TransAnno {
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
		// 情况1: snv包含了整个编码区，即整个编码区被删除
		// ...+++,,,+++...
		//   |---------|
		//|---------------|
		// 情况2: snv发生在整个编码区左边
		// |-|...+++,,,+++...
		//   |-|
		// 情况3: snv发生在整个编码区左边
		// ...+++,,,+++...|-|
		//              |-|
		if trans.Strand == "+" {
			transAnno.NAChange = fmt.Sprintf("c.%s_%sdelins%s", utrPosOfNAchange1, utrPosOfNAchange2, snv.Alt)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.%s_%sdelins%s", utrPosOfNAchange2, utrPosOfNAchange1, pkg.RevComp(snv.Alt))
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
			transAnno.NAChange = fmt.Sprintf("c.%s_%sdelins%s", utrPosOfNAchange1, cdsPosOfNAchange2, snv.Alt)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.%s_%sdelins%s", cdsPosOfNAchange2, utrPosOfNAchange1, pkg.RevComp(snv.Alt))
		}
	} else if trans.CdsEnd >= snv.Start && trans.CdsEnd < snv.End && trans.CdsStart <= snv.Start {
		// ...,,,+++...
		//        |--|
		//        |----|
		if trans.Strand == "+" {
			transAnno.NAChange = fmt.Sprintf("c.%s_%sdelins%s", cdsPosOfNAchange1, utrPosOfNAchange2, snv.Alt)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.%s_%sdelins%s", utrPosOfNAchange2, cdsPosOfNAchange1, pkg.RevComp(snv.Alt))
		}

	} else {
		// 变异跨region
		if !region1.Equal(region2) {
			if trans.Strand == "+" {
				transAnno.NAChange = fmt.Sprintf("c.%s_%sdelins%s", cdsPosOfNAchange1, cdsPosOfNAchange2, snv.Alt)
			} else {
				transAnno.NAChange = fmt.Sprintf("c.%s_%sdelins%s", cdsPosOfNAchange2, cdsPosOfNAchange1, pkg.RevComp(snv.Alt))
			}
		} else {
			if region1.Type == pkg.RType_INTRON {
				if trans.Strand == "+" {
					transAnno.NAChange = fmt.Sprintf("c.%s_%sdelins%s", cdsPosOfNAchange1, cdsPosOfNAchange2, snv.Alt)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%s_%sdelins%s", cdsPosOfNAchange2, cdsPosOfNAchange1, pkg.RevComp(snv.Alt))
				}
				if dist1 > 0 && dist2 > 0 && pkg.Min(dist1, dist2) <= 2 {
					transAnno.Event = "splicing"
					transAnno.Region = "splicing"
				}
			} else {
				cstart, cend := cLen1, cLen2
				cstart += snv.Start - region1.Start + 1
				cend += snv.End - region2.Start + 1
				transAnno = setSubAAChange(transAnno, trans, cstart, cend, snv.Alt)
			}
		}
	}
	return transAnno
}

func AnnoUnkSub(snv pkg.AnnoVariant, trans pkg.Transcript) TransAnno {
	alt := snv.Alt
	transAnno := NewTransAnno(trans)
	transAnno.Region2 = "ncRNA"
	nstart := snv.Start - trans.TxStart + 1
	nend := snv.End - trans.TxStart + 1
	dna := trans.DNA()
	ndna := pkg.Substitute2(dna, nstart, nend, alt)
	if trans.Strand == "-" {
		alt = pkg.RevComp(alt)
		dna = pkg.RevComp(dna)
		ndna = pkg.RevComp(ndna)
	}
	start := pkg.DifferenceSimple(dna, ndna)
	if nstart == nend {
		transAnno.NAChange = fmt.Sprintf("c.%ddelins%s", start, alt)
	} else {
		transAnno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", start, start+nend-nstart, alt)
	}
	return transAnno
}
