package genebased

import (
	"fmt"
	"open-anno/pkg/gene"
	"open-anno/pkg/seq"
	"open-anno/pkg/variant"
)

func findSnpRegion(regions gene.Regions, snv variant.Variant) (gene.Region, int) {
	var cLen int
	for _, region := range regions {
		if region.Start <= snv.Start && snv.End <= region.End {
			return region, cLen
		}
		if region.Type == gene.RType_CDS {
			cLen += region.End - region.Start + 1
		}
	}
	return gene.Region{}, cLen
}

func AnnoSnp(snv variant.Variant, trans gene.Transcript, aashort bool) SnvGeneBased {
	region, cLen := findSnpRegion(trans.Regions, snv)
	anno := NewSnvGeneBased(trans, region)
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
		if region.Type == gene.RType_INTRON {
			anno.Event = "."
			dist1 := snv.Start - region.Start + 1
			dist2 := region.End - snv.Start + 1
			if dist1 <= 2 || dist2 <= 2 {
				anno.Event = "splicing"
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
			start := seq.DifferenceSimple(cdna, ncdna)
			na1, na2 := cdna[start-1], ncdna[start-1]
			anno.NAChange = fmt.Sprintf("c.%d%c>%c", start, na1, na2)
			start = seq.DifferenceSimple(protein, protein)
			aa1, aa2 := protein[start-1], nprotein[start-1]
			if aa1 == aa2 {
				anno.Event = "synonymous_snv"
			} else {
				if aa1 == 'M' && start == 1 {
					anno.Event = "startloss"
				} else if aa1 == '*' {
					anno.Event = "stoploss"
				} else if aa2 == '*' {
					anno.Event = "stopgain"
				} else {
					anno.Event = "nonsynonymous_snv"
				}
			}
			anno.AAChange = fmt.Sprintf("p.%s%d%s", seq.AAName(aa1, aashort), start, seq.AAName(aa2, aashort))
		}
	}
	return anno
}
