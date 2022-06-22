package snv

import (
	"fmt"
	"open-anno/pkg/io"
	"open-anno/pkg/io/refgene"
	"open-anno/pkg/seq"
)

func findSnpRegion(regions refgene.Regions, snv io.Variant) (refgene.Region, int) {
	var cLen int
	for _, region := range regions {
		if region.Start <= snv.Start && snv.End <= region.End {
			return region, cLen
		}
		if region.Type == refgene.RType_CDS {
			cLen += region.End - region.Start + 1
		}
	}
	return refgene.Region{}, cLen
}

func AnnoSnp(snv io.Variant, trans refgene.Transcript) TransAnno {
	region, cLen := findSnpRegion(trans.Regions, snv)
	anno := NewTransAnno(trans, region)
	if region.End < trans.CdsStart {
		if trans.Strand == "+" {
			anno.NAChange = fmt.Sprintf("c.-%d%s>%s", trans.CdsStart-snv.End, snv.Ref, snv.Alt)
		} else {
			anno.NAChange = fmt.Sprintf("c.+%d%s>%s", trans.CdsStart-snv.End, seq.RevComp(snv.Ref), seq.RevComp(snv.Alt))
		}
	} else if region.Start > trans.CdsEnd {
		anno.Region = region.Name()
		if trans.Strand == "+" {
			anno.NAChange = fmt.Sprintf("c.+%d%s>%s", snv.Start-trans.CdsEnd, snv.Ref, snv.Alt)
		} else {
			anno.NAChange = fmt.Sprintf("c.-%d%s>%s", snv.Start-trans.CdsEnd, seq.RevComp(snv.Ref), seq.RevComp(snv.Alt))
		}
	} else {
		if region.Type == refgene.RType_INTRON {
			dist1 := snv.Start - region.Start + 1
			dist2 := region.End - snv.Start + 1
			if dist1 <= 2 || dist2 <= 2 {
				anno.Event = "splicing"
				anno.Region = "splicing"
			}
			if trans.Strand == "+" {
				if dist1 <= dist2 {
					anno.NAChange = fmt.Sprintf("c.%d+%d%s>%s", cLen, dist1, snv.Ref, snv.Alt)
				} else {
					anno.NAChange = fmt.Sprintf("c.%d-%d%s>%s", cLen+1, dist2, snv.Ref, snv.Alt)
				}
			} else {
				if dist1 <= dist2 {
					anno.NAChange = fmt.Sprintf("c.%d-%d%s>%s", trans.CLen()-cLen+1, dist1, seq.RevComp(snv.Ref), seq.RevComp(snv.Alt))
				} else {
					anno.NAChange = fmt.Sprintf("c.%d+%d%s>%s", trans.CLen()-cLen, dist2, seq.RevComp(snv.Ref), seq.RevComp(snv.Alt))
				}
			}

		} else {
			pos := cLen + snv.Start - region.Start + 1
			cdna := trans.CDNA()
			ncdna := seq.Substitute(cdna, pos, snv.Alt)
			if trans.Strand == "-" {
				cdna = seq.RevComp(cdna)
				ncdna = seq.RevComp(ncdna)
			}
			protein := seq.Translate(cdna, trans.Chrom == "MT")
			nprotein := seq.Translate(ncdna, trans.Chrom == "MT")
			cstart := seq.DifferenceSimple(cdna, ncdna)
			na1, na2 := cdna[cstart-1], ncdna[cstart-1]
			anno.NAChange = fmt.Sprintf("c.%d%c>%c", cstart, na1, na2)
			pstart := seq.DifferenceSimple(protein, protein)
			if pstart > len(protein) {
				pstart = cstart / 3

			}
			aa1, aa2 := protein[pstart-1], nprotein[pstart-1]
			if aa1 == aa2 {
				anno.Event = "synonymous_snv"
			} else {
				if aa1 == 'M' && pstart == 1 {
					anno.Event = "startloss"
				} else if aa1 == '*' {
					anno.Event = "stoploss"
				} else if aa2 == '*' {
					anno.Event = "stopgain"
				} else {
					anno.Event = "nonsynonymous_snv"
				}
			}
			anno.AAChange = fmt.Sprintf("p.%s%d%s", seq.AAName(aa1, AA_SHORT), pstart, seq.AAName(aa2, AA_SHORT))
		}
	}
	return anno
}
