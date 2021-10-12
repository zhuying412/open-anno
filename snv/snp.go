package snv

import (
	"fmt"
	"grandanno/gene"
	"grandanno/seq"
	"strings"
)

func GetCdsChangeChangeOfSNP(pos int, alt byte, cdna seq.Sequence, protein seq.Sequence, isMt bool) (event string, na_change string, aa_change string) {
	var varCdna, varProtein seq.Sequence
	varCdna = cdna.ChangeWithSnp(pos, alt)
	varProtein = varCdna.Translate(isMt)
	for i := 0; i < cdna.Len(); i++ {
		if j, na1, na2 := i/3, cdna.Base(i), varCdna.Base(i); na1 != na2 && j < protein.Len() {
			aa1, aa2 := protein.Base(j), varProtein.Base(j)
			if aa1 == aa2 {
				event = "synonymous_snp"
			} else {
				if aa1 == '*' {
					event = "stoploss"
				} else if aa2 == '*' {
					event = "stopgain"
				} else {
					event = "nonsynonymous_snp"
				}
			}
			na_change = fmt.Sprintf("c.%d%c>%c", i+1, na1, na2)
			aa_change = fmt.Sprintf("p.%s%d%s", seq.AAMap[aa1], j+1, seq.AAMap[aa2])
			break
		}
	}
	return event, na_change, aa_change
}

func NewAnnotationOfForwardSNP(snp Snv, refgene gene.Refgene) (anno Annotation) {
	cdna, protein := refgene.Cdna, refgene.Protein
	ref, alt := snp.Ref.Base(0), snp.Alt.Base(0)
	pos, regionCount := 0, len(refgene.Regions)
	for i := 0; i < regionCount; i++ {
		region := refgene.Regions[i]
		prevRegion, hasPrev := refgene.Regions.GetPrev(i, '+')
		nextRegion, hasNext := refgene.Regions.GetNext(i, '+')
		if region.Start > snp.Start {
			break
		} else if region.End < snp.Start {
			if region.Type == "cds" {
				pos += region.End - region.Start + 1
			}
		} else {
			if region.Type == "intron" {
				distance1 := snp.Start - region.Start + 1
				distance2 := region.End - snp.End + 1
				if distance1 < IntronDistance && hasPrev {
					if distance1 <= 2 {
						anno.Region = "splicing_site"
					} else {
						anno.Region = "splicing_region"
					}
					if prevRegion.Type == "cds" {
						if refgene.Tag == "cmpl" {
							anno.SetExon(prevRegion.ExonOrder)
							anno.NaChange = fmt.Sprintf("c.%d+%d%c>%c", pos, distance1, ref, alt)
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
							anno.NaChange = fmt.Sprintf("c.%d-%d%c>%c", pos+1, distance2, ref, alt)
						}
					} else {
						anno.Region = strings.Join([]string{nextRegion.Type, anno.Region}, "_")
					}
				} else {
					anno.Region = "intronic"
				}
			} else if strings.HasPrefix(region.Type, "utr") {
				if snp.Start-snp.Start < IntronDistance && hasPrev && prevRegion.Type == "intron" {
					anno.Region = strings.Join([]string{region.Type, "exon_splicing"}, "_")
				} else if region.End-snp.Start < IntronDistance && hasNext && nextRegion.Type == "intron" {
					anno.Region = strings.Join([]string{region.Type, "exon_splicing"}, "_")
				} else {
					anno.Region = region.Type
				}
			} else {
				if snp.Start-region.Start < IntronDistance && hasPrev && prevRegion.Type == "intron" {
					anno.Region = "CDS_splicing"
				} else if region.End-snp.Start < IntronDistance && hasNext && nextRegion.Type == "intron" {
					anno.Region = "CDS_splicing"
				} else {
					anno.Region = "exonic"
				}
				pos += snp.Start - region.Start + 1
				if refgene.Tag == "cmpl" {
					anno.SetExon(region.ExonOrder)
					anno.Event, anno.NaChange, anno.AaChange = GetCdsChangeChangeOfSNP(pos, alt, cdna, protein, refgene.Chrom == "MT")
				}
			}
		}
	}
	return anno
}

func NewAnnotationOfBackwardSNP(snp Snv, refgene gene.Refgene) (anno Annotation) {
	cdna, protein := refgene.Cdna, refgene.Protein
	ref, alt := snp.Ref.Base(0), snp.Alt.Base(0)
	pos, regionCount := 0, len(refgene.Regions)
	for i := 0; i < regionCount; i++ {
		region := refgene.Regions[i]
		prevRegion, hasPrev := refgene.Regions.GetPrev(i, '-')
		nextRegion, hasNext := refgene.Regions.GetNext(i, '-')
		if region.End < snp.Start {
			break
		} else if region.Start > snp.End {
			if region.Type == "cds" {
				pos += region.End - region.Start + 1
			}
		} else {
			if region.Type == "intron" {
				distance1 := snp.Start - region.Start + 1
				distance2 := region.End - snp.End + 1
				if distance1 <= IntronDistance && hasNext {
					if distance1 <= 2 {
						anno.Region = "splicing_site"
					} else {
						anno.Region = "splicing_region"
					}
					if nextRegion.Type == "cds" {
						if refgene.Tag == "cmpl" {
							anno.SetExon(nextRegion.ExonOrder)
							anno.NaChange = fmt.Sprintf("c.%d-%d%c>%c", pos+1, distance1, ref, alt)
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
							anno.NaChange = fmt.Sprintf("c.%d+%d%c>%c", pos, distance2, ref, alt)
						}
					} else {
						anno.Region = strings.Join([]string{prevRegion.Type, anno.Region}, "_")
					}
				} else {
					anno.Region = "intronic"
				}
			} else if strings.HasPrefix(region.Type, "utr") {
				if snp.Start-region.Start < IntronDistance && hasNext && nextRegion.Type == "intron" {
					anno.Region = strings.Join([]string{region.Type, "_exon_splicing"}, "_")
				} else if region.End-snp.Start < IntronDistance && hasPrev && prevRegion.Type == "intron" {
					anno.Region = strings.Join([]string{region.Type, "_exon_splicing"}, "_")
				} else {
					anno.Region = region.Type
				}
			} else {
				if snp.Start-region.Start < IntronDistance && hasNext && nextRegion.Type == "intron" {
					anno.Region = "CDS_splicing"
				} else if region.End-snp.Start < IntronDistance && hasPrev && prevRegion.Type == "intron" {
					anno.Region = "CDS_splicing"
				} else {
					anno.Region = "exonic"
				}
				pos += region.End - snp.Start + 1
				if refgene.Tag == "cmpl" {
					anno.SetExon(region.ExonOrder)
					anno.Event, anno.NaChange, anno.AaChange = GetCdsChangeChangeOfSNP(pos, alt, cdna, protein, refgene.Chrom == "MT")
				}
			}
		}
	}
	return anno
}

func NewAnnotationOfSNP(snp Snv, refgene gene.Refgene) Annotation {
	var anno Annotation
	if refgene.Strand == '+' {
		anno = NewAnnotationOfForwardSNP(snp, refgene)
	} else {
		anno = NewAnnotationOfBackwardSNP(snp, refgene)
	}
	anno.GeneSymbol = refgene.Gene
	anno.GeneEntrezId = refgene.EntrezId
	anno.Transcript = refgene.Transcript
	return anno
}
