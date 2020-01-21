package snv

import (
	"fmt"
	"grandanno/core"
	"strings"
)

func (anno *Annotation) annoCdsChangeOfSnp(pos int, alt byte, cdna core.Sequence, protein core.Sequence, isMt bool) {
	var varCdna, varProtein core.Sequence
	varCdna = cdna.GetSnpSequence(pos, alt)
	varProtein = varCdna.Translate(isMt)
	for i := 0; i < cdna.GetLen(); i++ {
		if j, na1, na2 := i/3, cdna.GetChar(i), varCdna.GetChar(i); na1 != na2 && j < protein.GetLen() {
			aa1, aa2 := protein.GetChar(j), varProtein.GetChar(j)
			if aa1 == aa2 {
				anno.Function = "synonymous_snv"
			} else {
				if aa1 == '*' {
					anno.Function = "stoploss"
				} else if aa2 == '*' {
					anno.Function = "stopgain"
				} else {
					anno.Function = "nonsynonymous_snv"
				}
			}
			anno.NaChange = fmt.Sprintf("c.%d%c>%c", i+1, na1, na2)
			anno.AaChange = fmt.Sprintf("p.%s%d%s", core.GetOne2Three(aa1), j+1, core.GetOne2Three(aa2))
			break
		}
	}
}

func (anno *Annotation) annoSnpForward(variant core.Variant, refgene core.Refgene, splicingLen int) {
	cdna, protein := refgene.Cdna, refgene.Protein
	ref, alt := variant.Ref.GetChar(0), variant.Alt.GetChar(0)
	pos, regionCount := 0, len(refgene.Regions)
	for i := 0; i < regionCount; i++ {
		region := refgene.Regions[i]
		prevRegion, hasPrev := refgene.Regions.GetPrev(i, '+')
		nextRegion, hasNext := refgene.Regions.GetNext(i, '+')
		if region.Start > variant.Start {
			break
		} else if region.End < variant.Start {
			if region.Typo == "cds" {
				pos += region.End - region.Start + 1
			}
		} else {
			if region.Typo == "intron" {
				distance1 := variant.Start - region.Start + 1
				distance2 := region.End - variant.End + 1
				if distance1 < splicingLen && hasPrev {
					if distance1 <= 2 {
						anno.Region = "splicing_site"
					} else {
						anno.Region = "splicing_region"
					}
					if prevRegion.Typo == "cds" {
						if refgene.Tag == "cmpl" {
							anno.SetExon(prevRegion.ExonOrder)
							anno.NaChange = fmt.Sprintf("c.%d+%d%c>%c", pos, distance1, ref, alt)
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
							anno.NaChange = fmt.Sprintf("c.%d-%d%c>%c", pos+1, distance2, ref, alt)
						}
					} else {
						anno.Region = strings.Join([]string{nextRegion.Typo, anno.Region}, "_")
					}
				} else {
					anno.Region = "intronic"
				}
			} else if strings.HasPrefix(region.Typo, "utr") {
				if variant.Start-variant.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = strings.Join([]string{region.Typo, "exon_splicing"}, "_")
				} else if region.End-variant.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = strings.Join([]string{region.Typo, "exon_splicing"}, "_")
				} else {
					anno.Region = region.Typo
				}
			} else {
				if variant.Start-region.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
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
						anno.annoCdsChangeOfSnp(pos, alt, cdna, protein, true)
					} else {
						anno.annoCdsChangeOfSnp(pos, alt, cdna, protein, false)
					}
				}
			}
		}
	}

}

func (anno Annotation) annoSnpBackward(variant core.Variant, refgene core.Refgene, splicingLen int) {
	cdna, protein := refgene.Cdna, refgene.Protein
	ref, alt := variant.Ref[0], variant.Alt[0]
	pos, regionCount := 0, len(refgene.Regions)
	for i := 0; i < regionCount; i++ {
		region := refgene.Regions[i]
		prevRegion, hasPrev := refgene.Regions.GetPrev(i, '-')
		nextRegion, hasNext := refgene.Regions.GetNext(i, '-')
		if region.End < variant.Start {
			break
		} else if region.Start > variant.End {
			if region.Typo == "cds" {
				pos += region.End - region.Start + 1
			}
		} else {
			if region.Typo == "intron" {
				distance1 := variant.Start - region.Start + 1
				distance2 := region.End - variant.End + 1
				if distance1 <= splicingLen && hasNext {
					if distance1 <= 2 {
						anno.Region = "splicing_site"
					} else {
						anno.Region = "splicing_region"
					}
					if nextRegion.Typo == "cds" {
						if refgene.Tag == "cmpl" {
							anno.SetExon(nextRegion.ExonOrder)
							anno.NaChange = fmt.Sprintf("c.%d-%d%c>%c", pos+1, distance1, ref, alt)
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
							anno.NaChange = fmt.Sprintf("c.%d+%d%c>%c", pos, distance2, ref, alt)
						}
					} else {
						anno.Region = strings.Join([]string{prevRegion.Typo, anno.Region}, "_")
					}
				} else {
					anno.Region = "intronic"
				}
			} else if strings.HasPrefix(region.Typo, "utr") {
				if variant.Start-region.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = strings.Join([]string{region.Typo, "_exon_splicing"}, "_")
				} else if region.End-variant.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = strings.Join([]string{region.Typo, "_exon_splicing"}, "_")
				} else {
					anno.Region = region.Typo
				}
			} else {
				if variant.Start-region.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else if region.End-variant.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else {
					anno.Region = "exonic"
				}
				pos += region.End - variant.Start + 1
				if refgene.Tag == "cmpl" {
					anno.SetExon(region.ExonOrder)
					anno.annoCdsChangeOfSnp(pos, alt, cdna, protein, refgene.Chrom == "MT")
				}
			}
		}
	}
}

func (anno *Annotation) AnnoSnp(snp Snv, refgene core.Refgene, splicingLen int) {
	if refgene.Strand == '+' {
		anno.annoSnpForward(snp.GetVariant(), refgene, splicingLen)
	} else {
		anno.annoSnpBackward(snp.GetVariant(), refgene, splicingLen)
	}
}
