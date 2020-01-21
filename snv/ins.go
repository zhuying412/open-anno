package snv

import (
	"fmt"
	"grandanno/core"
	"strings"
)

func (anno *Annotation) annoCdsChangeOfIns(pos int, alt core.Sequence, cdna core.Sequence, protein core.Sequence, isMt bool) {
	var varCdna, varProtein core.Sequence
	varCdna = cdna.GetInsSequence(pos, alt)
	if isMt {
		varProtein = cdna.Translate(true)
	} else {
		varProtein = cdna.Translate(false)
	}
	for i := 0; i < cdna.GetLen(); i++ {
		if cdna.GetChar(i) != varCdna.GetChar(i) {
			anno.NaChange = fmt.Sprintf("c.%dins%s", i, varCdna.GetSeq(i, len(alt)))
			break
		}
	}
	lenL, lenR := 0, 0
	for i := 0; i < protein.GetLen(); i++ {
		if protein.GetChar(i) != varProtein.GetChar(i) {
			break
		}
		lenL++
	}
	if lenL == protein.GetLen() {
		anno.Region = "utr3"
	} else {
		if alt.GetLen()%3 == 0 {
			for i := protein.GetLen() - 1; i >= 0; i-- {
				if protein.GetChar(i) != varProtein.GetChar(i) {
					break
				}
				lenR++
			}
			if lenL+lenR > protein.GetLen() {
				lenR = protein.GetLen() - lenL
			}
			start, end1, end2 := lenL, protein.GetLen()-lenR, varProtein.GetLen()-lenR
			altAa := varProtein.GetSeq(start, end2-start)
			if varProtein.GetIndex('*') < 0 {
				anno.Function = "ins_nonframeshift_stoploss"
			} else {
				if varProtein.GetIndex('*') < protein.GetLen()-1 {
					anno.Function = "ins_nonframeshift_stopgain"
				} else {
					anno.Function = "ins_nonframeshift"
				}
			}
			if start == end1 {
				anno.AaChange = fmt.Sprintf(
					"p.%s%d_%s%dins%s",
					core.GetOne2Three(protein.GetChar(start-1)),
					start,
					core.GetOne2Three(protein.GetChar(start)),
					start+1,
					altAa.GetOne2Tree(),
				)
			} else {
				refAa := protein.GetSeq(start, end1-start)
				anno.AaChange = fmt.Sprintf("p.%s%ddelins%s", refAa.GetOne2Tree(), start, altAa.GetOne2Tree())
			}
		} else {
			start := lenL
			if start < protein.GetLen() {
				if stopIndex := varProtein.GetIndex('*'); stopIndex > -1 && stopIndex < protein.GetLen()-1 {
					anno.Function = "ins_frameshift"
				} else {
					anno.Function = "ins_frameshift_stopgain"
				}
				anno.AaChange = fmt.Sprintf(
					"p.%s%d%sfs",
					core.GetOne2Three(protein.GetChar(start)),
					start+1,
					core.GetOne2Three(varProtein.GetChar(start)),
				)
			}
		}
	}
}

func (anno *Annotation) annoInsForward(variant core.Variant, refgene core.Refgene, splicingLen int) {
	cdna, protein := refgene.Cdna, refgene.Protein
	alt := variant.Alt
	pos, regionCount := 0, len(refgene.Regions)
	for i := 0; i < regionCount; i++ {
		region := refgene.Regions[i]
		prevRegion, hasPrev := refgene.Regions.GetPrev(i, '+')
		nextRegion, hasNext := refgene.Regions.GetNext(i, '+')
		if region.Start > variant.Start+1 {
			break
		} else if region.End <= variant.Start {
			if region.Typo == "cds" {
				pos += region.End - region.Start + 1
			}
		} else {
			if region.Typo == "intron" {
				distance1 := variant.Start - region.Start + 2
				distance2 := variant.Start - region.Start + 1
				if distance1 <= splicingLen && hasPrev {
					if distance1 <= 2 {
						anno.Region = "splicing_site"
					} else {
						anno.Region = "splicing_region"
					}
					if prevRegion.Typo == "cds" {
						if refgene.Tag == "cmpl" {
							anno.SetExon(prevRegion.ExonOrder)
							anno.NaChange = fmt.Sprintf("c.%d+%dins%s", pos, distance1, alt)
						}
					} else {
						anno.Region = strings.Join([]string{prevRegion.Typo, anno.Region}, "_")
					}
				} else if distance2 <= splicingLen && hasNext {
					if distance2 <= 2 {
						anno.Region = "splicing_site"
					} else {
						anno.Region = "splicing_region"
					}
					if nextRegion.Typo == "cds" {
						if refgene.Tag == "cmpl" {
							anno.SetExon(nextRegion.ExonOrder)
							anno.NaChange = fmt.Sprintf("c.%d-%dins%s", pos+1, distance2, alt)
						}
					} else {
						anno.Region = strings.Join([]string{nextRegion.Typo, anno.Region}, "_")
					}
				} else {
					anno.Region = "intronic"
				}
			} else if strings.HasPrefix(region.Typo, "utr") {
				if variant.Start-region.Start+1 < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = strings.Join([]string{region.Typo, "exon_splicing"}, "_")
				} else if region.End-variant.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = strings.Join([]string{region.Typo, "exon_splicing"}, "_")
				} else {
					anno.Region = region.Typo
				}
			} else {
				if variant.Start-region.Start+1 < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else if region.End-variant.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else {
					anno.Region = "exonic"
				}
				pos += variant.Start - region.Start + 1
				if refgene.Tag == "cmpl" {
					anno.SetExon(region.ExonOrder)
					if refgene.Chrom == "MT" {
						anno.annoCdsChangeOfIns(pos, alt, cdna, protein, true)
					} else {
						anno.annoCdsChangeOfIns(pos, alt, cdna, protein, false)
					}

				}
			}
		}
	}
}

func (anno *Annotation) annoInsBackward(variant core.Variant, refgene core.Refgene, splicingLen int) {
	cdna, protein := refgene.Cdna, refgene.Protein
	alt := variant.Alt
	pos, regionCount := 0, len(refgene.Regions)
	for i := regionCount - 1; i >= 0; i-- {
		region := refgene.Regions[i]
		prevRegion, hasPrev := refgene.Regions.GetPrev(i, '-')
		nextRegion, hasNext := refgene.Regions.GetNext(i, '-')
		if region.End < variant.Start {
			break
		} else if region.Start > variant.Start {
			if region.Typo == "cds" {
				pos += region.End - region.Start + 1
			}
		} else {
			if region.Typo == "intron" {
				distance1 := variant.Start - region.Start + 2
				distance2 := region.End - variant.Start + 1
				if distance1 <= splicingLen && hasNext {
					if distance1 <= 2 {
						anno.Region = "splicing_site"
					} else {
						anno.Region = "splicing_region"
					}
					if nextRegion.Typo == "cds" {
						if refgene.Tag == "cmpl" {
							anno.SetExon(nextRegion.ExonOrder)
							anno.NaChange = fmt.Sprintf("c.%d-%dins%s", pos+1, distance1, alt)
						}
					} else {
						anno.Region = strings.Join([]string{nextRegion.Typo, anno.Region}, "_")
					}
				} else if distance2 <= splicingLen && hasPrev {
					if distance2 <= 2 {
						anno.Region = "splicing_site"
					} else {
						anno.Region = "splicing_region"
					}
					if prevRegion.Typo == "cds" {
						if refgene.Tag == "cmpl" {
							anno.SetExon(prevRegion.ExonOrder)
							anno.NaChange = fmt.Sprintf("c.%d-%dins%s", pos, distance2, alt)
						}
					} else {
						anno.Region = strings.Join([]string{prevRegion.Typo, anno.Region}, "_")
					}
				} else {
					anno.Region = "intronic"
				}
			} else if strings.HasPrefix(region.Typo, "utr") {
				if variant.Start-region.Start+1 < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = strings.Join([]string{region.Typo, "exon_splicing"}, "_")
				} else if region.End-variant.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = strings.Join([]string{region.Typo, "exon_splicing"}, "_")
				} else {
					anno.Region = region.Typo
				}
			} else {
				if variant.Start-region.Start+1 < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else if region.End-variant.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else {
					anno.Region = "exonic"
				}
				pos += region.End - variant.Start
				if refgene.Tag == "cmpl" {
					anno.SetExon(region.ExonOrder)
					anno.annoCdsChangeOfIns(pos, alt, cdna, protein, refgene.Chrom == "MT")
				}
			}
		}
	}
}
func (anno *Annotation) AnnoIns(ins Snv, refgene core.Refgene, splicingLen int) {
	if refgene.Strand == '+' {
		anno.annoInsForward(ins.GetVariant(), refgene, splicingLen)
	} else {
		anno.annoInsBackward(ins.GetVariant(), refgene, splicingLen)
	}
}
