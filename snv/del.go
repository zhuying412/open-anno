package snv

import (
	"bytes"
	"grandanno/core"
	"strconv"
	"strings"
)

func (anno Annotation) setNaChangeOfSplicingDelOne(length int, flag byte, distance int) {
	var buffer bytes.Buffer
	buffer.WriteString("c.")
	buffer.WriteString(strconv.Itoa(length))
	buffer.WriteByte(flag)
	buffer.WriteString(strconv.Itoa(distance))
	buffer.WriteString("del")
	anno.NaChange = buffer.String()
}

func (anno Annotation) setNaChangeOfSplicingDelTwo(length int, flag byte, distance1 int, distance2 int) {
	var buffer bytes.Buffer
	buffer.WriteString("c.")
	buffer.WriteString(strconv.Itoa(length))
	buffer.WriteByte(flag)
	buffer.WriteString(strconv.Itoa(distance1))
	buffer.WriteByte('_')
	buffer.WriteString(strconv.Itoa(length))
	buffer.WriteByte(flag)
	buffer.WriteString(strconv.Itoa(distance2))
	buffer.WriteString("del")
	anno.NaChange = buffer.String()
}

func (anno Annotation) setNaChangeOfCdsDel(start int, end int) {
	var buffer bytes.Buffer
	buffer.WriteString("c.")
	buffer.WriteString(strconv.Itoa(start))
	buffer.WriteByte('_')
	buffer.WriteString(strconv.Itoa(end))
	buffer.WriteString("del")
	anno.NaChange = buffer.String()
}

func (anno Annotation) setAaChangeOfDelOne(ref byte, pos int) {
	var buffer bytes.Buffer
	buffer.WriteString("p.")
	buffer.WriteString(core.AaOne2ThreeDict[ref])
	buffer.WriteString(strconv.Itoa(pos))
	buffer.WriteString("del")
	anno.NaChange = buffer.String()
}

func (anno Annotation) setAaChangeOfDelMany(ref1 byte, start int, ref2 byte, end int) {
	var buffer bytes.Buffer
	buffer.WriteString("p.")
	buffer.WriteString(core.AaOne2ThreeDict[ref1])
	buffer.WriteString(strconv.Itoa(start))
	buffer.WriteByte('_')
	buffer.WriteString(core.AaOne2ThreeDict[ref2])
	buffer.WriteString(strconv.Itoa(end))
	buffer.WriteString("del")
	anno.NaChange = buffer.String()
}

func (anno Annotation) setAaChangeOfDelReplace(ref1 byte, start int, ref2 byte, end int, alt core.Sequence) {
	var buffer bytes.Buffer
	buffer.WriteString("p.")
	buffer.WriteString(core.AaOne2ThreeDict[ref1])
	buffer.WriteString(strconv.Itoa(start))
	buffer.WriteByte('_')
	buffer.WriteString(core.AaOne2ThreeDict[ref2])
	buffer.WriteString(strconv.Itoa(end))
	buffer.WriteString("insdel")
	for _, altChar := range alt {
		buffer.WriteString(core.AaOne2ThreeDict[altChar])
	}
	anno.NaChange = buffer.String()
}

func (anno Annotation) setAaChangeOfDelFrameshift(ref byte, pos int, alt byte) {
	var buffer bytes.Buffer
	buffer.WriteString("p.")
	buffer.WriteString(core.AaOne2ThreeDict[ref])
	buffer.WriteString(strconv.Itoa(pos))
	buffer.WriteString(core.AaOne2ThreeDict[alt])
	buffer.WriteString("fs")
	anno.NaChange = buffer.String()
}

func (anno Annotation) annoCdsChangeOfDel(lenL int, lenR int, cdna core.Sequence, protein core.Sequence, isMt bool) {
	var buffer bytes.Buffer
	if lenL > 0 {
		buffer.Write(cdna.GetSeq(0, lenL))
	}
	if lenR > 0 {
		buffer.Write(cdna.GetSeq(cdna.GetLen()-lenR, lenR))
	}
	var varCdna, varProtein core.Sequence
	varCdna = buffer.Bytes()
	varProtein = varCdna.Translate(isMt)
	lenl, lenr := 0, 0
	for i := 0; i < lenL+lenR; i++ {
		if cdna[i] != varCdna[i] {
			break
		}
		lenl++
	}
	for i := lenL + lenR - 1; i >= 0; i-- {
		if cdna[i] != varCdna[i] {
			break
		}
		lenr++
	}
	if lenl+lenr > lenL+lenR {
		lenr = lenL + lenR - lenl
	}
	anno.setNaChangeOfCdsDel(lenl+1, lenL+lenR-lenr)
	lenDel := lenL + lenR - lenl - lenr
	lenl, lenr, lenp, lenvp := 0, 0, protein.GetLen(), varProtein.GetLen()
	for i := 0; i < lenvp; i++ {
		if protein[i] == varProtein[i] {
			break
		}
		lenl++
	}
	if lenDel%3 == 0 {
		for i := lenvp - 1; i >= 0; i-- {
			if protein[i] != varProtein[i] {
				break
			}
			lenr++
		}
		if lenl+lenr > lenvp {
			lenr = lenvp - lenl
		}
		start, end1, end2 := lenl+1, lenp-lenr, lenvp-lenr
		if stopIndex := varProtein.GetIndex('*'); stopIndex < 0 {
			anno.Function = "del_nonframeshift_stoploss"
		} else if stopIndex < lenvp-1 {
			anno.Function = "del_nonframeshift_stopgain"
		} else {
			anno.Function = "del_nonframeshift"
		}
		if start == end2+1 {
			if start == end1 {
				anno.setAaChangeOfDelOne(protein.GetChar(start-1), start)
			} else {
				anno.setAaChangeOfDelMany(protein.GetChar(start-1), start, protein.GetChar(end1-1), end1)
			}
		} else {
			anno.setAaChangeOfDelReplace(protein.GetChar(start-1), start, protein.GetChar(end1-1), end1, varProtein.GetSeq(start-1, end2-start+1))
		}
	} else {
		start := lenl + 1
		if start > varProtein.GetLen() {
			anno.Function = "del_nonframeshift_stoploss"
			anno.setAaChangeOfDelMany(protein.GetChar(start-1), start, protein[lenp-1], lenp)
		} else {
			if stopIndex := varProtein.GetIndex('*'); stopIndex < 0 {
				anno.Function = "del_frameshift_stoploss"
			} else if stopIndex < lenvp-1 {
				anno.Function = "del_frameshift_stopgain"
			} else {
				anno.Function = "del_frameshift"
			}
			anno.setAaChangeOfDelFrameshift(protein.GetChar(start-1), start, varProtein.GetChar(start-1))
		}
	}
}

func (anno Annotation) annoIntronSplicingOfDel(del Snv, region core.Region, lenL int, ref core.Sequence, side byte, strand byte) {
	var distance1, distance2, length int
	var flag byte
	if side == 'l' {
		distance1 = del.Variant.Start - region.Start + 1
		distance2 = del.Variant.End - region.Start + 1
		if strand == '+' {
			length, flag = lenL, '+'
		} else {
			length, flag = lenL+1, '-'
		}
	} else {
		distance1 = region.End - del.Variant.Start + 1
		distance2 = region.End - del.Variant.End + 1
		if strand == '+' {
			length, flag = lenL+1, '-'
		} else {
			length, flag = lenL, '+'
		}
	}
	if ref.GetLen() == 1 {
		anno.setNaChangeOfSplicingDelOne(length, flag, distance1)
	} else {
		if flag == '+' && distance1 > distance2 || flag == '-' && distance1 < distance2 {
			distance1, distance2 = distance2, distance1
		}
		anno.setNaChangeOfSplicingDelTwo(length, flag, distance1, distance2)
	}
}

func (anno *Annotation) annoDelForward(del Snv, refgene core.Refgene, splicingLen int) {
	cdna, protein, ref := refgene.Cdna, refgene.Protein, del.Variant.Ref
	lenL, lenR := 0, 0
	regionCount := len(refgene.Regions)
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
		if region.Start > del.Variant.End {
			if region.Typo == "cds" {
				lenR += region.End - region.Start + 1
			}
		} else if region.End < del.Variant.Start {
			if region.Typo == "cds" {
				lenL += region.End - region.Start + 1
			}
		} else {
			if region.Start <= del.Variant.Start && del.Variant.End <= region.End {
				if region.Typo == "intron" {
					distance1 := del.Variant.Start - region.Start + 1
					distance2 := region.End - del.Variant.End + 1
					if distance1 <= splicingLen && hasPrev {
						if prevRegion.Typo == "cds" {
							if distance1 <= 2 {
								anno.Region = "splicing_site"
							} else {
								anno.Region = "splicing_region"
							}
							if refgene.Tag == "cmpl" {
								anno.SetExon(prevRegion.ExonOrder)
								anno.annoIntronSplicingOfDel(del, region, lenL, ref, 'l', '+')
							} else {
								if distance1 <= 2 {
									anno.Region = prevRegion.Typo + "_splicing_site"
								} else {
									anno.Region = prevRegion.Typo + "_splicing_region"
								}
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
								anno.annoIntronSplicingOfDel(del, region, lenL, ref, 'r', '+')
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
					if del.Variant.Start-region.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
						anno.Region = region.Typo + "_exon_splicing"
					} else if region.End-del.Variant.End < splicingLen && hasNext && nextRegion.Typo == "intron" {
						anno.Region = region.Typo + "_exon_splicing"
					} else {
						anno.Region = region.Typo
					}
				} else {
					lenL += del.Variant.Start - region.Start
					lenR += region.End - del.Variant.End
					if del.Variant.Start-region.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
						anno.Region = "CDS_splicing"
					} else if region.End-del.Variant.End < splicingLen && hasNext && nextRegion.Typo == "intron" {
						anno.Region = "CDS_splicing"
					} else {
						anno.Region = "exonic"
					}
				}
			} else {
				if region.Typo == "cds" {
					anno.Region = "exonic"
					anno.SetExon(region.ExonOrder)
					if region.Start < del.Variant.Start {
						lenL += del.Variant.Start - region.Start
						if hasNext && nextRegion.Typo == "intron" {
							anno.Region = "oCDS_splicing"
						}
					}
					if region.End > del.Variant.End {
						lenR += region.End - del.Variant.End
						if hasPrev && prevRegion.Typo == "intron" {
							anno.Region = "oCDS_splicing"
						}
					}
				} else {
					if anno.Region != "." {
						continue
					}
					if del.Variant.Start < region.Start && region.Start < del.Variant.End {
						if strings.HasPrefix(region.Typo, "utr") {
							if hasPrev {
								if prevRegion.Typo == "intron" {
									anno.Region = region.Typo + "_exon_splicing"
								} else {
									anno.Region = "exonic"
									anno.SetExon(prevRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Typo
							}
						} else {
							if hasPrev {
								if strings.HasPrefix(prevRegion.Typo, "utr") {
									anno.Region = prevRegion.Typo + "_exon_splicing"
								} else {
									anno.Region = "oCDS_splicing"
									anno.SetExon(prevRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Typo
							}
						}
					}
					if del.Variant.Start < region.End && region.End < del.Variant.End {
						if strings.HasPrefix(region.Typo, "utr") {
							if hasNext {
								if nextRegion.Typo == "intron" {
									anno.Region = region.Typo + "_exon_splicing"
								} else {
									anno.Region = "exonic"
									anno.SetExon(nextRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Typo
							}
						} else {
							if hasNext {
								if strings.HasPrefix(nextRegion.Typo, "utr") {
									anno.Region = nextRegion.Typo + "_exon_splicing"
								} else {
									anno.Region = "oCDS_splicing"
									anno.SetExon(nextRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Typo
							}
						}
					}
				}
			}
		}
	}
	if refgene.Tag == "cmpl" && (anno.Region == "exonic" || strings.HasSuffix(anno.Region, "CDS_splicing")) {
		anno.annoCdsChangeOfDel(lenL, lenR, cdna, protein, refgene.Chrom == "MT")
	}
}

func (anno *Annotation) annoDelBackward(del Snv, refgene core.Refgene, splicingLen int) {
	cdna, protein, ref := refgene.Cdna, refgene.Protein, del.Variant.Ref
	lenL, lenR := 0, 0
	regionCount := len(refgene.Regions)
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
		if region.End < del.Variant.Start {
			if region.Typo == "cds" {
				lenR += region.End - region.Start + 1
			}
		} else if region.Start > del.Variant.End {
			if region.Typo == "cds" {
				lenL += region.End - region.Start + 1
			}
		} else {
			if region.Start <= del.Variant.Start && del.Variant.End <= region.End {
				if region.Typo == "intron" {
					distance1 := del.Variant.Start - region.Start + 1
					distance2 := region.End - del.Variant.End + 1
					if distance1 <= splicingLen && hasNext {
						if nextRegion.Typo == "cds" {
							if distance1 <= 2 {
								anno.Region = "splicing_site"
							} else {
								anno.Region = "splicing_region"
							}
							if refgene.Tag == "cmpl" {
								anno.SetExon(nextRegion.ExonOrder)
								anno.annoIntronSplicingOfDel(del, region, lenL, ref, 'l', '-')
							} else {
								if distance1 <= 2 {
									anno.Region = nextRegion.Typo + "_splicing_site"
								} else {
									anno.Region = nextRegion.Typo + "_splicing_region"
								}
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
								anno.annoIntronSplicingOfDel(del, region, lenL, ref, 'r', '-')
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
					if del.Variant.Start-region.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
						anno.Region = region.Typo + "_exon_splicing"
					} else if region.End-del.Variant.End < splicingLen && hasNext && nextRegion.Typo == "intron" {
						anno.Region = region.Typo + "_exon_splicing"
					} else {
						anno.Region = region.Typo
					}
				} else {
					lenL += region.End - del.Variant.End
					lenR += del.Variant.Start - region.Start
					anno.SetExon(region.ExonOrder)
					if del.Variant.Start-region.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
						anno.Region = "CDS_splicing"
					} else if region.End-del.Variant.End < splicingLen && hasPrev && prevRegion.Typo == "intron" {
						anno.Region = "CDS_splicing"
					} else {
						anno.Region = "exonic"
					}
				}
			} else {
				if region.Typo == "cds" {
					anno.Region = "exonic"
					anno.SetExon(region.ExonOrder)
					if region.Start < del.Variant.Start {
						lenR += del.Variant.Start - region.Start
						if hasNext && nextRegion.Typo == "intron" {
							anno.Region = "oCDS_splicing"
						}
					}
					if region.End > del.Variant.End {
						lenL += region.End - del.Variant.End
						if hasPrev && prevRegion.Typo == "intron" {
							anno.Region = "oCDS_splicing"
						}
					}
				} else {
					if anno.Region != "." {
						continue
					}
					if del.Variant.Start < region.Start && region.Start < del.Variant.End {
						if strings.HasPrefix(region.Typo, "utr") {
							if hasNext {
								if nextRegion.Typo == "intron" {
									anno.Region = region.Typo + "_exon_splicing"
								} else {
									anno.Region = "exonic"
									anno.SetExon(nextRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Typo
							}
						} else {
							if hasNext {
								if strings.HasPrefix(nextRegion.Typo, "utr") {
									anno.Region = nextRegion.Typo + "_exon_splicing"
								} else {
									anno.Region = "oCDS_splicing"
									anno.SetExon(nextRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Typo
							}
						}
					}
					if del.Variant.Start < region.End && region.End < del.Variant.End {
						if strings.HasPrefix(region.Typo, "utr") {
							if strings.HasPrefix(region.Typo, "utr") {
								if hasPrev {
									if prevRegion.Typo == "intron" {
										anno.Region = region.Typo + "_exon_splicing"
									}
								} else {
									anno.Region = "exonic"
									anno.SetExon(prevRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Typo
							}
						} else {
							if strings.HasPrefix(region.Typo, "utr") {
								if hasPrev {
									if strings.HasPrefix(prevRegion.Typo, "utr") {
										anno.Region = prevRegion.Typo + "_exon_splicing"
									}
								} else {
									anno.Region = "oCDS_splicing"
									anno.SetExon(prevRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Typo
							}
						}
					}
				}
			}
		}
	}
	if refgene.Tag == "cmpl" && (anno.Region == "exonic" || strings.HasSuffix(anno.Region, "CDS_splicing")) {
		anno.annoCdsChangeOfDel(lenL, lenR, cdna, protein, refgene.Chrom == "MT")
	}
}
func (anno *Annotation) AnnoDel(del Snv, refgene core.Refgene, splicingLen int) {
	if refgene.Strand == '+' {
		anno.annoDelForward(del, refgene, splicingLen)
	} else {
		anno.annoDelBackward(del, refgene, splicingLen)
	}
}
