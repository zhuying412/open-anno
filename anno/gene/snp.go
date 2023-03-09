package gene

import (
	"fmt"
	"open-anno/anno"
	"open-anno/pkg"
)

func findSnpRegion(regions pkg.Regions, snv anno.AnnoVariant) (pkg.Region, int) {
	var cLen int
	for _, region := range regions {
		if region.Start <= snv.Start && snv.End <= region.End {
			return region, cLen
		}
		if region.Type == pkg.RType_CDS {
			cLen += region.End - region.Start + 1
		}
	}
	return pkg.Region{}, cLen
}

func AnnoSnp(snv anno.AnnoVariant, trans pkg.Transcript) TransAnno {
	region, cLen := findSnpRegion(trans.Regions, snv)
	transAnno := NewTransAnno(trans, region)
	if region.End < trans.CdsStart {
		if trans.Strand == "+" {
			transAnno.NAChange = fmt.Sprintf("c.-%d%s>%s", trans.CdsStart-snv.End, snv.Ref, snv.Alt)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.+%d%s>%s", trans.CdsStart-snv.End, pkg.RevComp(snv.Ref), pkg.RevComp(snv.Alt))
		}
	} else if region.Start > trans.CdsEnd {
		transAnno.Region = region.Name()
		if trans.Strand == "+" {
			transAnno.NAChange = fmt.Sprintf("c.+%d%s>%s", snv.Start-trans.CdsEnd, snv.Ref, snv.Alt)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.-%d%s>%s", snv.Start-trans.CdsEnd, pkg.RevComp(snv.Ref), pkg.RevComp(snv.Alt))
		}
	} else {
		if region.Type == pkg.RType_INTRON {
			dist1 := snv.Start - region.Start + 1
			dist2 := region.End - snv.Start + 1
			if dist1 <= 2 || dist2 <= 2 {
				transAnno.Event = "splicing"
				transAnno.Region = "splicing"
			}
			if trans.Strand == "+" {
				if dist1 <= dist2 {
					transAnno.NAChange = fmt.Sprintf("c.%d+%d%s>%s", cLen, dist1, snv.Ref, snv.Alt)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d-%d%s>%s", cLen+1, dist2, snv.Ref, snv.Alt)
				}
			} else {
				if dist1 <= dist2 {
					transAnno.NAChange = fmt.Sprintf("c.%d-%d%s>%s", trans.CLen()-cLen+1, dist1, pkg.RevComp(snv.Ref), pkg.RevComp(snv.Alt))
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d+%d%s>%s", trans.CLen()-cLen, dist2, pkg.RevComp(snv.Ref), pkg.RevComp(snv.Alt))
				}
			}

		} else {
			pos := cLen + snv.Start - region.Start + 1
			cdna := trans.CDNA()
			ncdna := pkg.Substitute(cdna, pos, snv.Alt)
			if trans.Strand == "-" {
				cdna = pkg.RevComp(cdna)
				ncdna = pkg.RevComp(ncdna)
			}
			protein := pkg.Translate(cdna, trans.Chrom == "MT")
			nprotein := pkg.Translate(ncdna, trans.Chrom == "MT")
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
