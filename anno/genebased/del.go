package genebased

import (
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/gene"
	"open-anno/pkg/seq"
	"open-anno/pkg/variant"
	"strings"
)

func findDelRegions(regions gene.Regions, snv variant.Variant) (gene.Regions, int) {
	var cLen int
	_regions := make(gene.Regions, 0)
	for _, region := range regions {
		if snv.Start <= region.End && snv.End >= region.Start {
			_regions = append(_regions, region)
		}
		if region.Type == gene.RType_CDS {
			if snv.End < region.Start {
				cLen += region.End - region.Start + 1
			}
		}
	}
	return _regions, cLen
}

func AnnoDel(snv variant.Variant, trans gene.Transcript, aashort bool) SnvGeneBased {
	var cLen int
	regions, cLen := findDelRegions(trans.Regions, snv)
	anno := NewSnvGeneBased(trans, regions[0])
	if snv.End < trans.CdsStart {
		if trans.Strand == "+" {
			anno.NAChange = fmt.Sprintf("c.-%d_-%ddel%s", trans.TxStart-snv.Start, trans.TxStart-snv.End, snv.Alt)
		} else {
			anno.NAChange = fmt.Sprintf("c.+%d_+%ddel%s", trans.TxStart-snv.End, trans.TxStart-snv.Start, seq.RevComp(snv.Alt))
		}
	} else if snv.Start > trans.CdsEnd {
		if trans.Strand == "+" {
			anno.NAChange = fmt.Sprintf("c.+%d_+%ddel%s", snv.Start-trans.TxEnd, snv.End-trans.TxEnd, snv.Alt)
		} else {
			anno.NAChange = fmt.Sprintf("c.-%d_-%ddel%s", snv.End-trans.TxEnd, snv.Start-trans.TxEnd, seq.RevComp(snv.Alt))
		}
	} else {
		hasNonCDS := false
		for _, r := range regions {
			if r.Type == gene.RType_INTRON {
				hasNonCDS = true
				break
			}
		}
		if hasNonCDS {
			if len(regions) == 1 {
				// 只有可能是Intron
				region := regions[0]
				anno.Event = "."
				dist1s, dist1e := snv.Start-region.Start+1, snv.End-region.Start+1
				dist2s, dist2e := region.End-region.Start+1, region.End-snv.End+1
				dist1, dist2 := dist1s, dist2e
				if pkg.Min(dist1, dist2) <= 2 {
					anno.Event = "splicing"
				}
				if trans.Strand == "+" {
					if dist1 <= dist2 {
						anno.NAChange = fmt.Sprintf("c.%d+%d_%d+%ddel%s", cLen, dist1s, cLen, dist1e, snv.Alt)
					} else {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%d-%ddel%s", cLen+1, dist2s, cLen+1, dist2e, snv.Alt)
					}
				} else {
					if dist1 <= dist2 {
						nclen := trans.CLen() - cLen + 1
						anno.NAChange = fmt.Sprintf("c.%d-%d_%d-%ddel%s", nclen, dist1e, nclen, dist1s, seq.RevComp(snv.Alt))
					} else {
						nclen := trans.CLen() - cLen
						anno.NAChange = fmt.Sprintf("c.%d+%d_%d+%ddel%s", nclen, dist1e, nclen, dist1s, seq.RevComp(snv.Alt))
					}
				}
			} else {
				// 1. start < txStart <= end <= txEnd
				// 2. txStart <= start <=  txEnd < end
				// 3. txStart <= start < end <= txEnd
			}
		} else {
			region := regions[0]
			start := cLen + snv.Start - region.Start + 1
			end := start + len(snv.Ref) - 1
			cdna := trans.CDNA()
			ncdna := seq.Delete(cdna, start, end)
			if trans.Strand == "-" {
				cdna = seq.RevComp(cdna)
				ncdna = seq.RevComp(ncdna)
			}
			protein := seq.Translate(cdna, trans.Chrom == "MT")
			nprotein := seq.Translate(ncdna, trans.Chrom == "MT")
			start = seq.DifferenceSimple(cdna, ncdna)
			alt := cdna[start-1 : start+len(snv.Alt)-1]
			anno.NAChange = fmt.Sprintf("c.%d_%ddel%s", start, start+len(snv.Alt)-1, alt)
			start, end1, end2 := seq.Difference(protein, nprotein)
			aa1 := protein[start-1 : end1]
			aa2 := nprotein[start-1 : end2]
			if len(snv.Alt)%3 == 0 {
				anno.Event = "del_nonframeshift"
				if len(aa2) == 0 {
					if len(aa1) == 1 {
						anno.AAChange = fmt.Sprintf("p.%s%ddel", seq.AAName(aa1, aashort), start)
					} else {
						anno.AAChange = fmt.Sprintf(
							"p.%s%d_%s%ddel",
							seq.AAName(aa1[0], aashort),
							start,
							seq.AAName(aa1[len(aa1)-1], aashort),
							end1,
						)
					}
				} else {
					if len(aa1) == 1 {
						anno.AAChange = fmt.Sprintf("p.%s%ddelins%s", seq.AAName(aa1, aashort), start, seq.AAName(aa2, aashort))
					} else {
						anno.AAChange = fmt.Sprintf(
							"p.%s%d_%s%ddelins%s",
							seq.AAName(aa1[0], aashort),
							start,
							seq.AAName(aa1[len(aa1)-1], aashort),
							end1,
							seq.AAName(aa2, aashort))
					}
				}
			} else {
				if start < len(protein) {
					anno.Event = "del_frameshift"
					var fs string
					fsi := strings.IndexByte(nprotein[start-1:], '*')
					if fsi == -1 {
						fs = "?"
					}
					if fsi == 0 {
						fs = fmt.Sprintf("%d", fsi)
					}
					anno.AAChange = fmt.Sprintf("p.%s%d%sfs*%s", seq.AAName(aa1[0], aashort), start, seq.AAName(aa2[0], aashort), fs)
				}
			}
		}
	}
	return anno
}
