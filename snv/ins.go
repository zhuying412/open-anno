package snv

import (
	"fmt"
	"grandanno/gene"
	"grandanno/seq"
	"strings"
)

func GetCdsChangeOfInsertion(pos int, alt seq.Sequence, cdna seq.Sequence, protein seq.Sequence, isMt bool) (event string, na_change string, aa_change string, region string) {
	var varCdna, varProtein seq.Sequence
	varCdna = cdna.ChangeWithIns(pos, alt)
	if isMt {
		varProtein = cdna.Translate(true)
	} else {
		varProtein = cdna.Translate(false)
	}
	for i := 0; i < cdna.Len(); i++ {
		if cdna.Base(i) != varCdna.Base(i) {
			aa_change = fmt.Sprintf("c.%dins%s", i, varCdna.SubSeq(i, len(alt)))
			break
		}
	}
	lenL, lenR := 0, 0
	for i := 0; i < protein.Len(); i++ {
		if protein.Base(i) != varProtein.Base(i) {
			break
		}
		lenL++
	}
	if lenL == protein.Len() {
		region = "utr3"
	} else {
		if alt.Len()%3 == 0 {
			for i := protein.Len() - 1; i >= 0; i-- {
				if protein.Base(i) != varProtein.Base(i) {
					break
				}
				lenR++
			}
			if lenL+lenR > protein.Len() {
				lenR = protein.Len() - lenL
			}
			start, end1, end2 := lenL, protein.Len()-lenR, varProtein.Len()-lenR
			altAa := varProtein.SubSeq(start, end2-start)
			if varProtein.Find('*') < 0 {
				event = "ins_nonframeshift_stoploss"
			} else {
				if varProtein.Find('*') < protein.Len()-1 {
					event = "ins_nonframeshift_stopgain"
				} else {
					event = "ins_nonframeshift"
				}
			}
			if start == end1 {
				aa_change = fmt.Sprintf(
					"p.%s%d_%s%dins%s",
					seq.AAMap[protein.Base(start-1)],
					start,
					seq.AAMap[protein.Base(start)],
					start+1,
					altAa.ProteinOne2Tree(),
				)
			} else {
				refAa := protein.SubSeq(start, end1-start)
				aa_change = fmt.Sprintf("p.%s%ddelins%s", refAa.ProteinOne2Tree(), start, altAa.ProteinOne2Tree())
			}
		} else {
			start := lenL
			if start < protein.Len() {
				if stopIndex := varProtein.Find('*'); stopIndex > -1 && stopIndex < protein.Len()-1 {
					event = "ins_frameshift"
				} else {
					event = "ins_frameshift_stopgain"
				}
				aa_change = fmt.Sprintf(
					"p.%s%d%sfs",
					seq.AAMap[protein.Base(start)],
					start+1,
					seq.AAMap[varProtein.Base(start)],
				)
			}
		}
	}
	return event, na_change, aa_change, region
}

func NewAnnotationOfForwardInsertion(ins Snv, refgene gene.Refgene) (anno Annotation) {
	cdna, protein := refgene.Cdna, refgene.Protein
	alt := ins.Alt
	pos, regionCount := 0, len(refgene.Regions)
	for i := 0; i < regionCount; i++ {
		region := refgene.Regions[i]
		prevRegion, hasPrev := refgene.Regions.GetPrev(i, '+')
		nextRegion, hasNext := refgene.Regions.GetNext(i, '+')
		if region.Start > ins.Start+1 {
			break
		} else if region.End <= ins.Start {
			if region.Type == "cds" {
				pos += region.End - region.Start + 1
			}
		} else {
			if region.Type == "intron" {
				distance1 := ins.Start - region.Start + 2
				distance2 := ins.Start - region.Start + 1
				if distance1 <= IntronDistance && hasPrev {
					if distance1 <= 2 {
						anno.Region = "splicing_site"
					} else {
						anno.Region = "splicing_region"
					}
					if prevRegion.Type == "cds" {
						if refgene.Tag == "cmpl" {
							anno.SetExon(prevRegion.ExonOrder)
							anno.NaChange = fmt.Sprintf("c.%d+%dins%s", pos, distance1, alt)
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
							anno.NaChange = fmt.Sprintf("c.%d-%dins%s", pos+1, distance2, alt)
						}
					} else {
						anno.Region = strings.Join([]string{nextRegion.Type, anno.Region}, "_")
					}
				} else {
					anno.Region = "intronic"
				}
			} else if strings.HasPrefix(region.Type, "utr") {
				if ins.Start-region.Start+1 < IntronDistance && hasPrev && prevRegion.Type == "intron" {
					anno.Region = strings.Join([]string{region.Type, "exon_splicing"}, "_")
				} else if region.End-ins.Start < IntronDistance && hasNext && nextRegion.Type == "intron" {
					anno.Region = strings.Join([]string{region.Type, "exon_splicing"}, "_")
				} else {
					anno.Region = region.Type
				}
			} else {
				if ins.Start-region.Start+1 < IntronDistance && hasPrev && prevRegion.Type == "intron" {
					anno.Region = "CDS_splicing"
				} else if region.End-ins.Start < IntronDistance && hasNext && nextRegion.Type == "intron" {
					anno.Region = "CDS_splicing"
				} else {
					anno.Region = "exonic"
				}
				pos += ins.Start - region.Start + 1
				if refgene.Tag == "cmpl" {
					anno.SetExon(region.ExonOrder)
					var newRegion string
					anno.Event, anno.NaChange, anno.AaChange, newRegion = GetCdsChangeOfInsertion(pos, alt, cdna, protein, refgene.Chrom == "MT")
					if newRegion == "" {
						anno.Region = newRegion
					}
				}
			}
		}
	}
	return anno
}

func NewAnnotationOfBackwardInsertion(ins Snv, refgene gene.Refgene) (anno Annotation) {
	cdna, protein := refgene.Cdna, refgene.Protein
	alt := ins.Alt
	pos, regionCount := 0, len(refgene.Regions)
	for i := regionCount - 1; i >= 0; i-- {
		region := refgene.Regions[i]
		prevRegion, hasPrev := refgene.Regions.GetPrev(i, '-')
		nextRegion, hasNext := refgene.Regions.GetNext(i, '-')
		if region.End < ins.Start {
			break
		} else if region.Start > ins.Start {
			if region.Type == "cds" {
				pos += region.End - region.Start + 1
			}
		} else {
			if region.Type == "intron" {
				distance1 := ins.Start - region.Start + 2
				distance2 := region.End - ins.Start + 1
				if distance1 <= IntronDistance && hasNext {
					if distance1 <= 2 {
						anno.Region = "splicing_site"
					} else {
						anno.Region = "splicing_region"
					}
					if nextRegion.Type == "cds" {
						if refgene.Tag == "cmpl" {
							anno.SetExon(nextRegion.ExonOrder)
							anno.NaChange = fmt.Sprintf("c.%d-%dins%s", pos+1, distance1, alt)
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
							anno.NaChange = fmt.Sprintf("c.%d-%dins%s", pos, distance2, alt)
						}
					} else {
						anno.Region = strings.Join([]string{prevRegion.Type, anno.Region}, "_")
					}
				} else {
					anno.Region = "intronic"
				}
			} else if strings.HasPrefix(region.Type, "utr") {
				if ins.Start-region.Start+1 < IntronDistance && hasNext && nextRegion.Type == "intron" {
					anno.Region = strings.Join([]string{region.Type, "exon_splicing"}, "_")
				} else if region.End-ins.Start < IntronDistance && hasPrev && prevRegion.Type == "intron" {
					anno.Region = strings.Join([]string{region.Type, "exon_splicing"}, "_")
				} else {
					anno.Region = region.Type
				}
			} else {
				if ins.Start-region.Start+1 < IntronDistance && hasNext && nextRegion.Type == "intron" {
					anno.Region = "CDS_splicing"
				} else if region.End-ins.Start < IntronDistance && hasPrev && prevRegion.Type == "intron" {
					anno.Region = "CDS_splicing"
				} else {
					anno.Region = "exonic"
				}
				pos += region.End - ins.Start
				if refgene.Tag == "cmpl" {
					anno.SetExon(region.ExonOrder)
					var newRegion string
					anno.Event, anno.NaChange, anno.AaChange, newRegion = GetCdsChangeOfInsertion(pos, alt, cdna, protein, refgene.Chrom == "MT")
					if newRegion == "" {
						anno.Region = newRegion
					}
				}
			}
		}
	}
	return anno
}

func NewAnnotationOfInsertion(ins Snv, refgene gene.Refgene) Annotation {
	var anno Annotation
	if refgene.Strand == '+' {
		anno = NewAnnotationOfForwardInsertion(ins, refgene)
	} else {
		anno = NewAnnotationOfBackwardInsertion(ins, refgene)
	}
	anno.GeneSymbol = refgene.Gene
	anno.GeneEntrezId = refgene.EntrezId
	anno.Transcript = refgene.Transcript
	return anno
}
