package snv

import (
	"bytes"
	"grandanno/core"
	"strconv"
	"strings"
)

func (anno *Annotation) setNaChangeOfSplicingIns(pos int, distance int, alt core.Sequence) {
	var buffer bytes.Buffer
	buffer.WriteString("c.")
	buffer.WriteString(strconv.Itoa(pos))
	if distance > 0 {
		buffer.WriteByte('+')
	}
	buffer.WriteString(strconv.Itoa(distance))
	buffer.WriteString("ins")
	buffer.Write(alt)
	anno.NaChange = buffer.String()
}

func (anno *Annotation) setNaChangeOfIns(pos int, alt core.Sequence) {
	var buffer bytes.Buffer
	buffer.WriteString("c.")
	buffer.WriteString(strconv.Itoa(pos))
	buffer.WriteString("ins")
	buffer.Write(alt)
	anno.NaChange = buffer.String()
}

func (anno *Annotation) setAaChangeOfInsBetween(ref1 byte, pos1 int, ref2 byte, pos2 int, alt core.Sequence) {
	var buffer bytes.Buffer
	buffer.WriteString("p.")
	buffer.WriteString(core.AaOne2ThreeDict[ref1])
	buffer.WriteString(strconv.Itoa(pos1))
	buffer.WriteByte('_')
	buffer.WriteString(core.AaOne2ThreeDict[ref2])
	buffer.WriteString(strconv.Itoa(pos2))
	buffer.WriteString("ins")
	for _, altChar := range alt {
		buffer.WriteString(core.AaOne2ThreeDict[altChar])
	}
	anno.NaChange = buffer.String()
}

func (anno *Annotation) setAaChangeOfInsReplace(ref core.Sequence, pos int, alt core.Sequence) {
	var buffer bytes.Buffer
	buffer.WriteString("p.")
	for _, refChar := range ref {
		buffer.WriteString(core.AaOne2ThreeDict[refChar])
	}
	buffer.WriteString(strconv.Itoa(pos))
	buffer.WriteString("delins")
	for _, altChar := range alt {
		buffer.WriteString(core.AaOne2ThreeDict[altChar])
	}
	anno.NaChange = buffer.String()
}

func (anno *Annotation) setAaChangeOfInsFrameshift(ref byte, pos int, alt byte) {
	var buffer bytes.Buffer
	buffer.WriteString("p.")
	buffer.WriteString(core.AaOne2ThreeDict[ref])
	buffer.WriteString(strconv.Itoa(pos))
	buffer.WriteString(core.AaOne2ThreeDict[alt])
	buffer.WriteString("fs")
	anno.NaChange = buffer.String()
}

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
			anno.setNaChangeOfIns(i, varCdna.GetSeq(i, len(alt)))
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
				anno.setAaChangeOfInsBetween(protein.GetChar(start-1), start, protein.GetChar(start), start+1, altAa)
			} else {
				refAa := protein.GetSeq(start, end1-start)
				anno.setAaChangeOfInsReplace(refAa, start, altAa)
			}
		} else {
			start := lenL
			if start < protein.GetLen() {
				if stopIndex := varProtein.GetIndex('*'); stopIndex > -1 && stopIndex < protein.GetLen()-1 {
					anno.Function = "ins_frameshift"
				} else {
					anno.Function = "ins_frameshift_stopgain"
				}
				anno.setAaChangeOfInsFrameshift(protein.GetChar(start), start+1, varProtein.GetChar(start))
			}
		}
	}
}

func (anno *Annotation) annoInsForward(ins Snv, refgene core.Refgene, splicingLen int) {
	cdna, protein := refgene.Cdna, refgene.Protein
	alt := ins.Variant.Alt
	pos, regionCount := 0, len(refgene.Regions)
	for i := 0; i < regionCount; i++ {
		region := refgene.Regions[i]
		hasPrev, hasNext := false, false
		var prevRegion, nextRegion core.Region
		if i-1 >= 0 {
			hasPrev = true
			prevRegion = refgene.Regions[i-1]
		}
		if i+1 < regionCount {
			hasNext = true
			nextRegion = refgene.Regions[i+1]
		}
		if region.Start > ins.Variant.Start+1 {
			break
		} else if region.End <= ins.Variant.Start {
			if region.Typo == "cds" {
				pos += region.End - region.Start + 1
			}
		} else {
			if region.Typo == "intron" {
				distance1 := ins.Variant.Start - region.Start + 2
				distance2 := ins.Variant.Start - region.Start + 1
				if distance1 <= splicingLen && hasPrev {
					if prevRegion.Typo == "cds" {
						if distance1 <= 2 {
							anno.Region = "splicing_site"
						} else {
							anno.Region = "splicing_region"
						}
						if refgene.Tag == "cmpl" {
							anno.SetExon(prevRegion.ExonOrder)
							anno.setNaChangeOfSplicingIns(pos, distance1, alt)
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
							anno.setNaChangeOfSplicingIns(pos+1, -distance2, alt)
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
				if ins.Variant.Start-region.Start+1 < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = region.Typo + "_exon_splicing"
				} else if region.End-ins.Variant.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = region.Typo + "_exon_splicing"
				} else {
					anno.Region = region.Typo
				}
			} else {
				if ins.Variant.Start-region.Start+1 < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else if region.End-ins.Variant.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else {
					anno.Region = "exonic"
				}
				pos += ins.Variant.Start - region.Start + 1
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

func (anno *Annotation) annoInsBackward(ins Snv, refgene core.Refgene, splicingLen int) {
	cdna, protein := refgene.Cdna, refgene.Protein
	alt := ins.Variant.Alt
	pos, regionCount := 0, len(refgene.Regions)
	for i := regionCount - 1; i >= 0; i-- {
		region := refgene.Regions[i]
		hasPrev, hasNext := false, false
		var prevRegion, nextRegion core.Region
		if i-1 >= 0 {
			hasNext = true
			nextRegion = refgene.Regions[i-1]
		}
		if i+1 < regionCount {
			hasPrev = true
			prevRegion = refgene.Regions[i+1]
		}
		if region.End < ins.Variant.Start {
			break
		} else if region.Start > ins.Variant.Start {
			if region.Typo == "cds" {
				pos += region.End - region.Start + 1
			}
		} else {
			if region.Typo == "intron" {
				distance1 := ins.Variant.Start - region.Start + 2
				distance2 := region.End - ins.Variant.Start + 1
				if distance1 <= splicingLen && hasNext {
					if nextRegion.Typo == "cds" {
						if distance1 <= 2 {
							anno.Region = "splicing_site"
						} else {
							anno.Region = "splicing_region"
						}
						if refgene.Tag == "cmpl" {
							anno.SetExon(nextRegion.ExonOrder)
							anno.setNaChangeOfSplicingIns(pos+1, -distance1, alt)
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
							distance := region.End - ins.Variant.Start + 1
							anno.setNaChangeOfSplicingIns(pos, -distance, alt)
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
				if ins.Variant.Start-region.Start+1 < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = region.Typo + "_exon_splicing"
				} else if region.End-ins.Variant.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = region.Typo + "_exon_splicing"
				} else {
					anno.Region = region.Typo
				}
			} else {
				if ins.Variant.Start-region.Start+1 < splicingLen && hasNext && nextRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else if region.End-ins.Variant.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
					anno.Region = "CDS_splicing"
				} else {
					anno.Region = "exonic"
				}
				pos += region.End - ins.Variant.Start
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
func (anno *Annotation) AnnoIns(ins Snv, refgene core.Refgene, splicingLen int) {
	if refgene.Strand == '+' {
		anno.annoInsForward(ins, refgene, splicingLen)
	} else {
		anno.annoInsBackward(ins, refgene, splicingLen)
	}
}
