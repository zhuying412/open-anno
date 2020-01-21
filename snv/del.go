package snv

import (
	"fmt"
	"grandanno/core"
	"strings"
)

func (anno Annotation) annoCdsChangeOfDel(lenL int, lenR int, cdna core.Sequence, protein core.Sequence, isMt bool) {
	varCdna := cdna.GetDelSequence(lenL, lenR)
	varProtein := varCdna.Translate(isMt)
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
	anno.NaChange = fmt.Sprintf("c.%d_%ddel", lenl+1, lenL+lenR-lenr)
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
				anno.AaChange = fmt.Sprintf("p.%s%ddel", core.GetOne2Three(protein.GetChar(start-1)), start)
			} else {
				anno.AaChange = fmt.Sprintf(
					"p.%s%d_%s%ddel",
					core.GetOne2Three(protein.GetChar(start-1)),
					start,
					core.GetOne2Three(protein.GetChar(end1-1)),
					end1,
				)
			}
		} else {
			anno.AaChange = fmt.Sprintf(
				"p.%s%d_%s%dinsdel%s",
				core.GetOne2Three(protein.GetChar(start-1)),
				start,
				core.GetOne2Three(protein.GetChar(end1-1)),
				end1,
				varProtein.GetSeq(start-1, end2-start+1).GetOne2Tree(),
			)
		}
	} else {
		start := lenl + 1
		if start > varProtein.GetLen() {
			anno.Function = "del_nonframeshift_stoploss"
			anno.AaChange = fmt.Sprintf(
				"p.%s%d_%s%s%ddel",
				core.GetOne2Three(protein.GetChar(start-1)),
				start,
				core.GetOne2Three(protein[lenp-1]),
				lenp,
			)
		} else {
			if stopIndex := varProtein.GetIndex('*'); stopIndex < 0 {
				anno.Function = "del_frameshift_stoploss"
			} else if stopIndex < lenvp-1 {
				anno.Function = "del_frameshift_stopgain"
			} else {
				anno.Function = "del_frameshift"
			}
			anno.AaChange = fmt.Sprintf(
				"p.%s%d%sfs",
				core.GetOne2Three(protein.GetChar(start-1)),
				start,
				core.GetOne2Three(varProtein.GetChar(start-1)),
			)
		}
	}
}

func (anno Annotation) annoIntronSplicingOfDel(variant core.Variant, region core.Region, lenL int, ref core.Sequence, side byte, strand byte) {
	var distance1, distance2, length int
	var flag byte
	if side == 'l' {
		distance1 = variant.Start - region.Start + 1
		distance2 = variant.End - region.Start + 1
		if strand == '+' {
			length, flag = lenL, '+'
		} else {
			length, flag = lenL+1, '-'
		}
	} else {
		distance1 = region.End - variant.Start + 1
		distance2 = region.End - variant.End + 1
		if strand == '+' {
			length, flag = lenL+1, '-'
		} else {
			length, flag = lenL, '+'
		}
	}
	if ref.GetLen() == 1 {
		anno.NaChange = fmt.Sprintf("c.%d%c%ddel", length, flag, distance1)
	} else {
		if flag == '+' && distance1 > distance2 || flag == '-' && distance1 < distance2 {
			distance1, distance2 = distance2, distance1
		}
		anno.NaChange = fmt.Sprintf("c.%d%c%d_%d%c%ddel", length, flag, distance1, length, flag, distance2)
	}
}

func (anno *Annotation) annoDelForward(variant core.Variant, refgene core.Refgene, splicingLen int) {
	cdna, protein, ref := refgene.Cdna, refgene.Protein, variant.Ref
	lenL, lenR := 0, 0
	regionCount := len(refgene.Regions)
	for i := 0; i < regionCount; i++ {
		region := refgene.Regions[i]
		prevRegion, hasPrev := refgene.Regions.GetPrev(i, '+')
		nextRegion, hasNext := refgene.Regions.GetNext(i, '+')
		if region.Start > variant.End {
			if region.Typo == "cds" {
				lenR += region.End - region.Start + 1
			}
		} else if region.End < variant.Start {
			if region.Typo == "cds" {
				lenL += region.End - region.Start + 1
			}
		} else {
			if region.Start <= variant.Start && variant.End <= region.End {
				if region.Typo == "intron" {
					distance1 := variant.Start - region.Start + 1
					distance2 := region.End - variant.End + 1
					if distance1 <= splicingLen && hasPrev {
						if distance1 <= 2 {
							anno.Region = "splicing_site"
						} else {
							anno.Region = "splicing_region"
						}
						if prevRegion.Typo == "cds" {
							if refgene.Tag == "cmpl" {
								anno.SetExon(prevRegion.ExonOrder)
								anno.annoIntronSplicingOfDel(variant, region, lenL, ref, 'l', '+')
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
								anno.annoIntronSplicingOfDel(variant, region, lenL, ref, 'r', '+')
							}
						} else {
							anno.Region = strings.Join([]string{nextRegion.Typo, anno.Region}, "_")
						}
					} else {
						anno.Region = "intronic"
					}
				} else if strings.HasPrefix(region.Typo, "utr") {
					if variant.Start-region.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
						anno.Region = strings.Join([]string{region.Typo, "exon_splicing"}, "_")
					} else if region.End-variant.End < splicingLen && hasNext && nextRegion.Typo == "intron" {
						anno.Region = strings.Join([]string{region.Typo, "exon_splicing"}, "_")
					} else {
						anno.Region = region.Typo
					}
				} else {
					lenL += variant.Start - region.Start
					lenR += region.End - variant.End
					if variant.Start-region.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
						anno.Region = "CDS_splicing"
					} else if region.End-variant.End < splicingLen && hasNext && nextRegion.Typo == "intron" {
						anno.Region = "CDS_splicing"
					} else {
						anno.Region = "exonic"
					}
				}
			} else {
				if region.Typo == "cds" {
					anno.Region = "exonic"
					anno.SetExon(region.ExonOrder)
					if region.Start < variant.Start {
						lenL += variant.Start - region.Start
						if hasNext && nextRegion.Typo == "intron" {
							anno.Region = "oCDS_splicing"
						}
					}
					if region.End > variant.End {
						lenR += region.End - variant.End
						if hasPrev && prevRegion.Typo == "intron" {
							anno.Region = "oCDS_splicing"
						}
					}
				} else {
					if anno.Region != "." {
						continue
					}
					if variant.Start < region.Start && region.Start < variant.End {
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
					if variant.Start < region.End && region.End < variant.End {
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

func (anno *Annotation) annoDelBackward(variant core.Variant, refgene core.Refgene, splicingLen int) {
	cdna, protein, ref := refgene.Cdna, refgene.Protein, variant.Ref
	lenL, lenR := 0, 0
	regionCount := len(refgene.Regions)
	for i := regionCount - 1; i >= 0; i-- {
		region := refgene.Regions[i]
		prevRegion, hasPrev := refgene.Regions.GetPrev(i, '+')
		nextRegion, hasNext := refgene.Regions.GetNext(i, '+')
		if region.End < variant.Start {
			if region.Typo == "cds" {
				lenR += region.End - region.Start + 1
			}
		} else if region.Start > variant.End {
			if region.Typo == "cds" {
				lenL += region.End - region.Start + 1
			}
		} else {
			if region.Start <= variant.Start && variant.End <= region.End {
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
								anno.annoIntronSplicingOfDel(variant, region, lenL, ref, 'l', '-')
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
								anno.annoIntronSplicingOfDel(variant, region, lenL, ref, 'r', '-')
							}
						} else {
							anno.Region = strings.Join([]string{prevRegion.Typo, anno.Region}, "_")
						}
					} else {
						anno.Region = "intronic"
					}
				} else if strings.HasPrefix(region.Typo, "utr") {
					if variant.Start-region.Start < splicingLen && hasPrev && prevRegion.Typo == "intron" {
						anno.Region = strings.Join([]string{region.Typo, "exon_splicing"}, "_")
					} else if region.End-variant.End < splicingLen && hasNext && nextRegion.Typo == "intron" {
						anno.Region = strings.Join([]string{region.Typo, "exon_splicing"}, "_")
					} else {
						anno.Region = region.Typo
					}
				} else {
					lenL += region.End - variant.End
					lenR += variant.Start - region.Start
					anno.SetExon(region.ExonOrder)
					if variant.Start-region.Start < splicingLen && hasNext && nextRegion.Typo == "intron" {
						anno.Region = "CDS_splicing"
					} else if region.End-variant.End < splicingLen && hasPrev && prevRegion.Typo == "intron" {
						anno.Region = "CDS_splicing"
					} else {
						anno.Region = "exonic"
					}
				}
			} else {
				if region.Typo == "cds" {
					anno.Region = "exonic"
					anno.SetExon(region.ExonOrder)
					if region.Start < variant.Start {
						lenR += variant.Start - region.Start
						if hasNext && nextRegion.Typo == "intron" {
							anno.Region = "oCDS_splicing"
						}
					}
					if region.End > variant.End {
						lenL += region.End - variant.End
						if hasPrev && prevRegion.Typo == "intron" {
							anno.Region = "oCDS_splicing"
						}
					}
				} else {
					if anno.Region != "." {
						continue
					}
					if variant.Start < region.Start && region.Start < variant.End {
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
					if variant.Start < region.End && region.End < variant.End {
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
		anno.annoDelForward(del.GetVariant(), refgene, splicingLen)
	} else {
		anno.annoDelBackward(del.GetVariant(), refgene, splicingLen)
	}
}
