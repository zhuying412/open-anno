package snv

import (
	"fmt"
	"grandanno/gene"
	"grandanno/seq"
	"strings"
)

func GetCdsChangeOfDeletion(lenL int, lenR int, cdna seq.Sequence, protein seq.Sequence, isMt bool) (event string, na_change string, aa_change string) {
	varCdna := cdna.ChangeWithDel(lenL, lenR)
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
	na_change = fmt.Sprintf("c.%d_%ddel", lenl+1, lenL+lenR-lenr)
	lenDel := lenL + lenR - lenl - lenr
	lenl, lenr, lenp, lenvp := 0, 0, protein.Len(), varProtein.Len()
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
		if stopIndex := varProtein.Find('*'); stopIndex < 0 {
			event = "del_nonframeshift_stoploss"
		} else if stopIndex < lenvp-1 {
			event = "del_nonframeshift_stopgain"
		} else {
			event = "del_nonframeshift"
		}
		if start == end2+1 {
			if start == end1 {
				aa_change = fmt.Sprintf("p.%s%ddel", seq.AAMap[protein.Base(start-1)], start)
			} else {
				aa_change = fmt.Sprintf(
					"p.%s%d_%s%ddel",
					seq.AAMap[protein.Base(start-1)],
					start,
					seq.AAMap[protein.Base(end1-1)],
					end1,
				)
			}
		} else {
			aa_change = fmt.Sprintf(
				"p.%s%d_%s%dinsdel%s",
				seq.AAMap[protein.Base(start-1)],
				start,
				seq.AAMap[protein.Base(end1-1)],
				end1,
				varProtein.SubSeq(start-1, end2-start+1).ProteinOne2Tree(),
			)
		}
	} else {
		start := lenl + 1
		if start > varProtein.Len() {
			aa_change = "del_nonframeshift_stoploss"
			aa_change = fmt.Sprintf(
				"p.%s%d_%s%s%ddel",
				seq.AAMap[protein.Base(start-1)],
				start,
				seq.AAMap[protein[lenp-1]],
				lenp,
			)
		} else {
			if stopIndex := varProtein.Find('*'); stopIndex < 0 {
				event = "del_frameshift_stoploss"
			} else if stopIndex < lenvp-1 {
				event = "del_frameshift_stopgain"
			} else {
				event = "del_frameshift"
			}
			event = fmt.Sprintf(
				"p.%s%d%sfs",
				seq.AAMap[protein.Base(start-1)],
				start,
				seq.AAMap[varProtein.Base(start-1)],
			)
		}
	}
	return event, na_change, aa_change
}

func GetIntronSplicingOfDeletion(del Snv, region gene.Region, lenL int, ref seq.Sequence, side byte, strand byte) (na_change string) {
	var distance1, distance2, length int
	var flag byte
	if side == 'l' {
		distance1 = del.Start - region.Start + 1
		distance2 = del.End - region.Start + 1
		if strand == '+' {
			length, flag = lenL, '+'
		} else {
			length, flag = lenL+1, '-'
		}
	} else {
		distance1 = region.End - del.Start + 1
		distance2 = region.End - del.End + 1
		if strand == '+' {
			length, flag = lenL+1, '-'
		} else {
			length, flag = lenL, '+'
		}
	}
	if ref.Len() == 1 {
		na_change = fmt.Sprintf("c.%d%c%ddel", length, flag, distance1)
	} else {
		if flag == '+' && distance1 > distance2 || flag == '-' && distance1 < distance2 {
			distance1, distance2 = distance2, distance1
		}
		na_change = fmt.Sprintf("c.%d%c%d_%d%c%ddel", length, flag, distance1, length, flag, distance2)
	}
	return na_change
}

func NewAnnotationOfForwardDeletion(del Snv, refgene gene.Refgene) (anno Annotation) {
	cdna, protein, ref := refgene.Cdna, refgene.Protein, del.Ref
	lenL, lenR := 0, 0
	regionCount := len(refgene.Regions)
	for i := 0; i < regionCount; i++ {
		region := refgene.Regions[i]
		prevRegion, hasPrev := refgene.Regions.GetPrev(i, '+')
		nextRegion, hasNext := refgene.Regions.GetNext(i, '+')
		if region.Start > del.End {
			if region.Type == "cds" {
				lenR += region.End - region.Start + 1
			}
		} else if region.End < del.Start {
			if region.Type == "cds" {
				lenL += region.End - region.Start + 1
			}
		} else {
			if region.Start <= del.Start && del.End <= region.End {
				if region.Type == "intron" {
					distance1 := del.Start - region.Start + 1
					distance2 := region.End - del.End + 1
					if distance1 <= IntronDistance && hasPrev {
						if distance1 <= 2 {
							anno.Region = "splicing_site"
						} else {
							anno.Region = "splicing_region"
						}
						if prevRegion.Type == "cds" {
							if refgene.Tag == "cmpl" {
								anno.SetExon(prevRegion.ExonOrder)
								anno.NaChange = GetIntronSplicingOfDeletion(del, region, lenL, ref, 'l', '+')
							}
						} else {
							anno.Region = strings.Join([]string{prevRegion.Type, anno.Region}, "_")
						}
					} else if distance2 <= IntronDistance && hasNext {
						if distance2 <= 2 {
							anno.Region = "splicing_site"
						} else {
							anno.Region = "splicing_region"
						}
						if nextRegion.Type == "cds" {
							if refgene.Tag == "cmpl" {
								anno.SetExon(nextRegion.ExonOrder)
								anno.NaChange = GetIntronSplicingOfDeletion(del, region, lenL, ref, 'r', '+')
							}
						} else {
							anno.Region = strings.Join([]string{nextRegion.Type, anno.Region}, "_")
						}
					} else {
						anno.Region = "intronic"
					}
				} else if strings.HasPrefix(region.Type, "utr") {
					if del.Start-region.Start < IntronDistance && hasPrev && prevRegion.Type == "intron" {
						anno.Region = strings.Join([]string{region.Type, "exon_splicing"}, "_")
					} else if region.End-del.End < IntronDistance && hasNext && nextRegion.Type == "intron" {
						anno.Region = strings.Join([]string{region.Type, "exon_splicing"}, "_")
					} else {
						anno.Region = region.Type
					}
				} else {
					lenL += del.Start - region.Start
					lenR += region.End - del.End
					if del.Start-region.Start < IntronDistance && hasPrev && prevRegion.Type == "intron" {
						anno.Region = "CDS_splicing"
					} else if region.End-del.End < IntronDistance && hasNext && nextRegion.Type == "intron" {
						anno.Region = "CDS_splicing"
					} else {
						anno.Region = "exonic"
					}
				}
			} else {
				if region.Type == "cds" {
					anno.Region = "exonic"
					anno.SetExon(region.ExonOrder)
					if region.Start < del.Start {
						lenL += del.Start - region.Start
						if hasNext && nextRegion.Type == "intron" {
							anno.Region = "oCDS_splicing"
						}
					}
					if region.End > del.End {
						lenR += region.End - del.End
						if hasPrev && prevRegion.Type == "intron" {
							anno.Region = "oCDS_splicing"
						}
					}
				} else {
					if anno.Region != "." {
						continue
					}
					if del.Start < region.Start && region.Start < del.End {
						if strings.HasPrefix(region.Type, "utr") {
							if hasPrev {
								if prevRegion.Type == "intron" {
									anno.Region = region.Type + "_exon_splicing"
								} else {
									anno.Region = "exonic"
									anno.SetExon(prevRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Type
							}
						} else {
							if hasPrev {
								if strings.HasPrefix(prevRegion.Type, "utr") {
									anno.Region = prevRegion.Type + "_exon_splicing"
								} else {
									anno.Region = "oCDS_splicing"
									anno.SetExon(prevRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Type
							}
						}
					}
					if del.Start < region.End && region.End < del.End {
						if strings.HasPrefix(region.Type, "utr") {
							if hasNext {
								if nextRegion.Type == "intron" {
									anno.Region = region.Type + "_exon_splicing"
								} else {
									anno.Region = "exonic"
									anno.SetExon(nextRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Type
							}
						} else {
							if hasNext {
								if strings.HasPrefix(nextRegion.Type, "utr") {
									anno.Region = nextRegion.Type + "_exon_splicing"
								} else {
									anno.Region = "oCDS_splicing"
									anno.SetExon(nextRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Type
							}
						}
					}
				}
			}
		}
	}
	if refgene.Tag == "cmpl" && (anno.Region == "exonic" || strings.HasSuffix(anno.Region, "CDS_splicing")) {
		anno.Event, anno.NaChange, anno.AaChange = GetCdsChangeOfDeletion(lenL, lenR, cdna, protein, refgene.Chrom == "MT")
	}
	return anno
}

func NewAnnotationOfBackwardDeletion(del Snv, refgene gene.Refgene) (anno Annotation) {
	cdna, protein, ref := refgene.Cdna, refgene.Protein, del.Ref
	lenL, lenR := 0, 0
	regionCount := len(refgene.Regions)
	for i := regionCount - 1; i >= 0; i-- {
		region := refgene.Regions[i]
		prevRegion, hasPrev := refgene.Regions.GetPrev(i, '+')
		nextRegion, hasNext := refgene.Regions.GetNext(i, '+')
		if region.End < del.Start {
			if region.Type == "cds" {
				lenR += region.End - region.Start + 1
			}
		} else if region.Start > del.End {
			if region.Type == "cds" {
				lenL += region.End - region.Start + 1
			}
		} else {
			if region.Start <= del.Start && del.End <= region.End {
				if region.Type == "intron" {
					distance1 := del.Start - region.Start + 1
					distance2 := region.End - del.End + 1
					if distance1 <= IntronDistance && hasNext {
						if distance1 <= 2 {
							anno.Region = "splicing_site"
						} else {
							anno.Region = "splicing_region"
						}
						if nextRegion.Type == "cds" {
							if refgene.Tag == "cmpl" {
								anno.SetExon(nextRegion.ExonOrder)
								anno.NaChange = GetIntronSplicingOfDeletion(del, region, lenL, ref, 'l', '-')
							}
						} else {
							anno.Region = strings.Join([]string{nextRegion.Type, anno.Region}, "_")
						}
					} else if distance2 <= IntronDistance && hasPrev {
						if distance2 <= 2 {
							anno.Region = "splicing_site"
						} else {
							anno.Region = "splicing_region"
						}
						if prevRegion.Type == "cds" {
							if refgene.Tag == "cmpl" {
								anno.SetExon(prevRegion.ExonOrder)
								anno.NaChange = GetIntronSplicingOfDeletion(del, region, lenL, ref, 'r', '-')
							}
						} else {
							anno.Region = strings.Join([]string{prevRegion.Type, anno.Region}, "_")
						}
					} else {
						anno.Region = "intronic"
					}
				} else if strings.HasPrefix(region.Type, "utr") {
					if del.Start-region.Start < IntronDistance && hasPrev && prevRegion.Type == "intron" {
						anno.Region = strings.Join([]string{region.Type, "exon_splicing"}, "_")
					} else if region.End-del.End < IntronDistance && hasNext && nextRegion.Type == "intron" {
						anno.Region = strings.Join([]string{region.Type, "exon_splicing"}, "_")
					} else {
						anno.Region = region.Type
					}
				} else {
					lenL += region.End - del.End
					lenR += del.Start - region.Start
					anno.SetExon(region.ExonOrder)
					if del.Start-region.Start < IntronDistance && hasNext && nextRegion.Type == "intron" {
						anno.Region = "CDS_splicing"
					} else if region.End-del.End < IntronDistance && hasPrev && prevRegion.Type == "intron" {
						anno.Region = "CDS_splicing"
					} else {
						anno.Region = "exonic"
					}
				}
			} else {
				if region.Type == "cds" {
					anno.Region = "exonic"
					anno.SetExon(region.ExonOrder)
					if region.Start < del.Start {
						lenR += del.Start - region.Start
						if hasNext && nextRegion.Type == "intron" {
							anno.Region = "oCDS_splicing"
						}
					}
					if region.End > del.End {
						lenL += region.End - del.End
						if hasPrev && prevRegion.Type == "intron" {
							anno.Region = "oCDS_splicing"
						}
					}
				} else {
					if anno.Region != "." {
						continue
					}
					if del.Start < region.Start && region.Start < del.End {
						if strings.HasPrefix(region.Type, "utr") {
							if hasNext {
								if nextRegion.Type == "intron" {
									anno.Region = region.Type + "_exon_splicing"
								} else {
									anno.Region = "exonic"
									anno.SetExon(nextRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Type
							}
						} else {
							if hasNext {
								if strings.HasPrefix(nextRegion.Type, "utr") {
									anno.Region = nextRegion.Type + "_exon_splicing"
								} else {
									anno.Region = "oCDS_splicing"
									anno.SetExon(nextRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Type
							}
						}
					}
					if del.Start < region.End && region.End < del.End {
						if strings.HasPrefix(region.Type, "utr") {
							if strings.HasPrefix(region.Type, "utr") {
								if hasPrev {
									if prevRegion.Type == "intron" {
										anno.Region = region.Type + "_exon_splicing"
									}
								} else {
									anno.Region = "exonic"
									anno.SetExon(prevRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Type
							}
						} else {
							if strings.HasPrefix(region.Type, "utr") {
								if hasPrev {
									if strings.HasPrefix(prevRegion.Type, "utr") {
										anno.Region = prevRegion.Type + "_exon_splicing"
									}
								} else {
									anno.Region = "oCDS_splicing"
									anno.SetExon(prevRegion.ExonOrder)
								}
							} else {
								anno.Region = region.Type
							}
						}
					}
				}
			}
		}
	}
	if refgene.Tag == "cmpl" && (anno.Region == "exonic" || strings.HasSuffix(anno.Region, "CDS_splicing")) {
		anno.Event, anno.NaChange, anno.AaChange = GetCdsChangeOfDeletion(lenL, lenR, cdna, protein, refgene.Chrom == "MT")
	}
	return anno
}

func NewAnnotationOfDeletion(del Snv, refgene gene.Refgene) Annotation {
	var anno Annotation
	if refgene.Strand == '+' {
		anno = NewAnnotationOfForwardDeletion(del, refgene)
	} else {
		anno = NewAnnotationOfBackwardDeletion(del, refgene)
	}
	anno.GeneSymbol = refgene.Gene
	anno.GeneEntrezId = refgene.EntrezId
	anno.Transcript = refgene.Transcript
	return anno
}
