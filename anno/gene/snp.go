package gene

import (
	"fmt"
	"open-anno/pkg"
)

func getUTRPosOfNAchange(trans pkg.Transcript, utrLen1, utrLen2, pos, uLen int, region pkg.Region) string {
	if pos < trans.CdsStart {
		if trans.Strand == "+" {
			if utrLen1 == 0 {
				return fmt.Sprintf("-%d", trans.CdsStart-pos)
			} else {
				if region.Type == pkg.RType_UTR {
					return fmt.Sprintf("-%d", utrLen1-(uLen+pos-region.Start))
				} else if region.Type == pkg.RType_INTRON {
					dist1, dist2 := pos-region.Start+1, region.End-pos+1
					if dist1 < dist2 {
						return fmt.Sprintf("-%d+%d", utrLen1-uLen+1, dist1)
					} else {
						return fmt.Sprintf("-%d-%d", utrLen1-uLen, dist2)
					}
				} else {
					return fmt.Sprintf("-%d-%d", utrLen1, trans.TxStart-pos)
				}
			}
		} else {
			if utrLen1 == 0 {
				return fmt.Sprintf("*%d", trans.CdsStart-pos)
			} else {
				if region.Type == pkg.RType_UTR {
					return fmt.Sprintf("*%d", utrLen1-(uLen+pos-region.Start))
				} else if region.Type == pkg.RType_INTRON {
					dist1, dist2 := pos-region.Start+1, region.End-pos+1
					if dist1 < dist2 {
						return fmt.Sprintf("*%d-%d", utrLen1-uLen+1, dist1)
					} else {
						return fmt.Sprintf("*%d+%d", utrLen1-uLen, dist2)
					}
				} else {
					return fmt.Sprintf("*%d+%d", utrLen1, trans.TxStart-pos)
				}
			}
		}
	}
	if pos > trans.CdsEnd {
		if trans.Strand == "+" {
			if utrLen2 == 0 {
				return fmt.Sprintf("*%d", pos-trans.CdsEnd)
			} else {
				if region.Type == pkg.RType_UTR {
					return fmt.Sprintf("*%d", uLen+(pos-region.Start+1))
				} else if region.Type == pkg.RType_INTRON {
					dist1, dist2 := pos-region.Start+1, region.End-pos+1
					if dist1 < dist2 {
						return fmt.Sprintf("*%d+%d", uLen, dist1)
					} else {
						return fmt.Sprintf("*%d-%d", uLen+1, dist2)
					}
				} else {
					return fmt.Sprintf("*%d+%d", utrLen2, pos-trans.TxEnd)
				}
			}
		} else {
			if utrLen2 == 0 {
				return fmt.Sprintf("-%d", pos-trans.CdsEnd)
			} else {
				if region.Type == pkg.RType_UTR {
					return fmt.Sprintf("-%d", uLen+(pos-region.Start+1))
				} else if region.Type == pkg.RType_INTRON {
					dist1, dist2 := pos-region.Start+1, region.End-pos+1
					if dist1 < dist2 {
						return fmt.Sprintf("-%d-%d", uLen, dist1)
					} else {
						return fmt.Sprintf("-%d+%d", uLen+1, dist2)
					}
				} else {
					return fmt.Sprintf("-%d-%d", utrLen2, pos-trans.TxEnd)
				}
			}
		}
	}
	return ""
}

func AnnoSnp(snv pkg.AnnoVariant, trans pkg.Transcript) TransAnno {
	region, cLen, uLen := trans.Region(snv.Start)
	transAnno := NewTransAnno(trans, region)
	if snv.Start < trans.CdsStart || snv.Start > trans.CdsEnd {
		utrLen1, utrLen2 := trans.ULen()
		utrPosOfNAchange := getUTRPosOfNAchange(trans, utrLen1, utrLen2, snv.Start, uLen, region)
		if trans.Strand == "+" {
			transAnno.NAChange = fmt.Sprintf("c.%s%s>%s", utrPosOfNAchange, snv.Ref, snv.Alt)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.%s%s>%s", utrPosOfNAchange, pkg.RevComp(snv.Ref), pkg.RevComp(snv.Alt))
		}
	} else {
		cdsLen := trans.CLen()
		if region.Type == pkg.RType_INTRON {
			dist1, dist2 := snv.Start-region.Start+1, region.End-snv.Start+1
			if trans.Strand == "+" {

				if dist1 < dist2 {
					transAnno.NAChange = fmt.Sprintf("c.%d+%d%s>%s", cLen, dist1, snv.Ref, snv.Alt)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d-%d%s>%s", cLen+1, dist2, snv.Ref, snv.Alt)
				}
			} else {
				if dist1 < dist2 {
					transAnno.NAChange = fmt.Sprintf("c.%d-%d%s>%s", cdsLen-cLen+1, dist1, pkg.RevComp(snv.Ref), pkg.RevComp(snv.Alt))
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d+%d%s>%s", cdsLen-cLen, dist2, pkg.RevComp(snv.Ref), pkg.RevComp(snv.Alt))
				}
			}
			if pkg.Min(dist1, dist2) <= 2 {
				transAnno.Event = "splicing"
				transAnno.Region = "splicing"
			}
		} else {
			pos := cLen + snv.Start - region.Start + 1
			cdna := trans.CDNA()
			ncdna := pkg.Substitute(cdna, pos, snv.Alt)
			if trans.Strand == "-" {
				cdna = pkg.RevComp(cdna)
				ncdna = pkg.RevComp(ncdna)
			}
			protein := pkg.Translate(cdna, trans.Chrom == "MT" || trans.Chrom == "chrM")
			nprotein := pkg.Translate(ncdna, trans.Chrom == "MT" || trans.Chrom == "chrM")
			cstart := pkg.DifferenceSimple(cdna, ncdna)
			na1, na2 := cdna[cstart-1], ncdna[cstart-1]
			transAnno.NAChange = fmt.Sprintf("c.%d%c>%c", cstart, na1, na2)
			pstart := pkg.DifferenceSimple(protein, nprotein)
			if pstart > len(protein) {
				pstart = cstart / 3

			}
			aa1, aa2 := protein[pstart-1], nprotein[pstart-1]
			if aa1 == aa2 {
				transAnno.Event = "synonymous"
			} else {
				if aa1 == 'M' && pstart == 1 {
					transAnno.Event = "startloss"
				} else if aa1 == '*' {
					transAnno.Event = "stoploss"
				} else if aa2 == '*' {
					transAnno.Event = "nonsense"
				} else {
					transAnno.Event = "missense"
				}
			}
			if aa1 == '*' {
				transAnno.AAChange = fmt.Sprintf("p.%s%d%sext*?", pkg.AAName(aa1, AA_SHORT), pstart, pkg.AAName(aa2, AA_SHORT))
			} else {
				transAnno.AAChange = fmt.Sprintf("p.%s%d%s", pkg.AAName(aa1, AA_SHORT), pstart, pkg.AAName(aa2, AA_SHORT))
			}

		}
	}
	return transAnno
}

func AnnoUnkSnp(snv pkg.AnnoVariant, trans pkg.Transcript) TransAnno {
	transAnno := NewTransAnno(trans)
	transAnno.Region = "ncRNA"
	pos := snv.Start - trans.TxStart + 1
	dna := trans.DNA()
	ndna := pkg.Substitute(dna, pos, snv.Alt)
	if trans.Strand == "-" {
		dna = pkg.RevComp(dna)
		ndna = pkg.RevComp(ndna)
	}
	start := pkg.DifferenceSimple(dna, ndna)
	na1, na2 := dna[start-1], ndna[start-1]
	transAnno.NAChange = fmt.Sprintf("n.%d%c>%c", start, na1, na2)
	return transAnno
}
