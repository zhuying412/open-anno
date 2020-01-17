package snv

import (
	"bytes"
	"grandanno/core"
	"strconv"
	"strings"
)

func (anno *Annotation) setNaChangeOfCdsSnp(pos int, refChar byte, altChar byte) {
	var buffer bytes.Buffer
	buffer.WriteString("c.")
	buffer.WriteString(strconv.Itoa(pos))
	buffer.Write([]byte{refChar, '>', altChar})
	anno.NaChange = buffer.String()
}

func (anno *Annotation) setNaChangeOfSplicingSnp(pos int, distance int, refChar byte, altChar byte) {
	var buffer bytes.Buffer
	buffer.WriteString("c.")
	buffer.WriteString(strconv.Itoa(pos))
	if distance > 0 {
		buffer.WriteByte('+')
	}
	buffer.WriteString(strconv.Itoa(distance))
	buffer.Write([]byte{refChar, '>', altChar})
	anno.NaChange = buffer.String()
}

func (anno *Annotation) setAaChangeOfCdsSnp(pos int, refChar byte, altChar byte) {
	var buffer bytes.Buffer
	buffer.WriteString("p.")
	buffer.WriteString(core.AaOne2ThreeDict[refChar])
	buffer.WriteString(strconv.Itoa(pos))
	buffer.WriteString(core.AaOne2ThreeDict[altChar])
	anno.AaChange = buffer.String()
}

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
			anno.setNaChangeOfCdsSnp(i+1, na1, na2)
			anno.setAaChangeOfCdsSnp(j+1, aa1, aa2)
			break
		}
	}
}

func (anno *Annotation) annoSnpForward(snp Snv, refgene core.Refgene, splicingLen int) {
	cdna, protein := refgene.Cdna, refgene.Protein
	ref, alt := snp.Variant.Ref.GetChar(0), snp.Variant.Alt.GetChar(0)
	pos, regionCount := 0, len(refgene.Regions)
	for i := 0; i < regionCount; i++ {
		region := refgene.Regions[i]
		hasPrev, prevRegion, hasNext, nextRegion := refgene.Regions.GetPreNextRegion(i, true)
		if region.Start > snp.Variant.Start {
			break
		} else if region.End < snp.Variant.Start {
			if region.Typo == "cds" {
				pos += region.End - region.Start + 1
			}
		} else {
			if region.Typo == "intron" {
				distance1 := snp.Variant.Start - region.Start + 1
				distance2 := region.End - snp.Variant.Start + 1
				if distance1 < splicingLen && hasPrev {
					if prevRegion.Typo == "cds" {
						if distance1 <= 2 {
							anno.Region = "splicing_site"
						} else {
							anno.Region = "splicing_region"
						}
						if refgene.Tag == "cmpl" {
							anno.SetExon(prevRegion.ExonOrder)
							anno.setNaChangeOfSplicingSnp(pos, distance1, ref, alt)
						}
					} else {
						if distance1 <= 2 {
							anno.Region = prevRegion.Typo + "_splicing_site"
						} else {
							anno.Region = prevRegion.Typo + "_splicing_region"
						}
					}
				} else if distance2 <= splicingLen && hasNext {
					if nextRegion.Typo == "cds" {
						if distance2 <= 2 {
							anno.Region = "splicing_site"
						} else {
							anno.Region = "splicing_region"
						}
						if refgene.Tag == "cmpl" {
							anno.SetExon(nextRegion.ExonOrder)
							anno.setNaChangeOfSplicingSnp(pos+1, -distance2, ref, alt)
						}
					} else {
						if distance2 <= 2 {
							anno.Region = nextRegion.Typo + "_splicing_site"
						} else {
							anno.Region = nextRegion.Typo + "_splicing_region"
						}
					}
				} else {
					anno.Region = "intronic"
				}
			} else if strings.HasPrefix(region.Typo, "utr") {
				if snp.Variant.Start-snp.Variant.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = region.Typo + "_exon_splicing"
				} else if region.End-snp.Variant.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = region.Typo + "_exon_splicing"
				} else {
					anno.Region = region.Typo
				}
			} else {
				if snp.Variant.Start-region.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else if region.End-snp.Variant.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else {
					anno.Region = "exonic"
				}
				pos += snp.Variant.Start - region.Start + 1
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

func (anno Annotation) annoSnpBackward(snp Snv, refgene core.Refgene, splicingLen int) {
	cdna, protein := refgene.Cdna, refgene.Protein
	ref, alt := snp.Variant.Ref[0], snp.Variant.Alt[0]
	pos, regionCount := 0, len(refgene.Regions)
	for i := 0; i < regionCount; i++ {
		region := refgene.Regions[i]
		hasPrev, prevRegion, hasNext, nextRegion := refgene.Regions.GetPreNextRegion(i, false)
		if region.End < snp.Variant.Start {
			break
		} else if region.Start > snp.Variant.Start {
			if region.Typo == "cds" {
				pos += region.End - region.Start + 1
			}
		} else {
			if region.Typo == "intron" {
				distance1 := snp.Variant.Start - region.Start + 1
				distance2 := region.End - snp.Variant.Start + 1
				if distance1 <= splicingLen && hasNext {
					if nextRegion.Typo == "cds" {
						if distance1 <= 2 {
							anno.Region = "splicing_site"
						} else {
							anno.Region = "splicing_region"
						}
						if refgene.Tag == "cmpl" {
							anno.SetExon(nextRegion.ExonOrder)
							anno.setNaChangeOfSplicingSnp(pos+1, -distance1, ref, alt)
						}
					} else {
						if distance1 <= 2 {
							anno.Region = nextRegion.Typo + "_splicing_site"
						} else {
							anno.Region = nextRegion.Typo + "_splicing_region"
						}
					}
				} else if distance2 <= splicingLen && hasPrev {
					if prevRegion.Typo == "cds" {
						if distance2 <= 2 {
							anno.Region = "splicing_site"
						} else {
							anno.Region = "splicing_region"
						}
						if refgene.Tag == "cmpl" {
							anno.SetExon(prevRegion.ExonOrder)
							anno.setNaChangeOfSplicingSnp(pos, distance2, ref, alt)
						}
					} else {
						if distance2 <= 2 {
							anno.Region = prevRegion.Typo + "_splicing_site"
						} else {
							anno.Region = prevRegion.Typo + "_splicing_region"
						}
					}
				} else {
					anno.Region = "intronic"
				}
			} else if strings.HasPrefix(region.Typo, "utr") {
				if snp.Variant.Start-region.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = region.Typo + "_exon_splicing"
				} else if region.End-snp.Variant.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = region.Typo + "_exon_splicing"
				} else {
					anno.Region = region.Typo
				}
			} else {
				if snp.Variant.Start-region.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else if region.End-snp.Variant.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else {
					anno.Region = "exonic"
				}
				pos += region.End - snp.Variant.Start + 1
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
		anno.annoSnpForward(snp, refgene, splicingLen)
	} else {
		anno.annoSnpBackward(snp, refgene, splicingLen)
	}
}
