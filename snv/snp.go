package snv

import (
	"fmt"
	"grandanno/gene"
	"grandanno/seq"
	"strings"
)

func NewAnnotationOfIntronSnp(snp Snv, refgene gene.Refgene, regionIndex int, exonLen int) (anno Annotation) {
	region := refgene.Regions[regionIndex]
	prevRegion, _ := refgene.Regions.GetPrev(regionIndex, refgene.Strand)
	nextRegion, _ := refgene.Regions.GetNext(regionIndex, refgene.Strand)
	distance1 := snp.Start - region.Start + 1
	distance2 := region.End - snp.End + 1
	var closestRegion gene.Region
	var pos, distance int
	var flag byte
	if distance1 <= distance2 {
		distance = distance1
	} else {
		distance = distance2
	}
	if (refgene.Strand == '+') == (distance1 <= distance2) {
		closestRegion = prevRegion
		pos = exonLen
		flag = '+'
	} else {
		closestRegion = nextRegion
		pos = exonLen + 1
		flag = '-'
	}
	anno.Region = "intronic"
	if distance <= SplicingDistance {
		anno.Region = "splicing"
		if !closestRegion.IsCDS() {
			anno.Region = closestRegion.Type + "_splicing"
		}
	}
	if closestRegion.IsCDS() {
		anno.SetExon(closestRegion.ExonOrder)
		if snp.Type() == "snp" {
			anno.NaChange = fmt.Sprintf("c.%d%c%d%s>%s", pos, flag, distance, snp.Ref, snp.Alt)
		} else {
			anno.NaChange = fmt.Sprintf("c.%d%c%dins%s", pos, flag, distance, snp.Alt)
		}
	}
	return anno
}

func NewAnnotationOfUtrSnp(refgene gene.Refgene, regionIndex int) (anno Annotation) {
	region := refgene.Regions[regionIndex]
	anno.Region = region.Type
	return Annotation{Region: region.Type}
}

func NewAnnotationOfCdsSnp(snp Snv, refgene gene.Refgene, regionIndex int, exonLen int) (anno Annotation) {
	region := refgene.Regions[regionIndex]
	var pos int
	anno.SetExon(region.ExonOrder)
	anno.Region = "exonic"
	if refgene.Strand == '+' {
		pos = exonLen + snp.Start - region.Start + 1
	} else {
		pos = exonLen + region.End - snp.Start + 1
	}
	if refgene.IsCmpl() {
		cdna, protein := refgene.Cdna, refgene.Protein
		newCdna := cdna.ChangeWithSnp(pos, snp.Alt.Base(0))
		newProtein := newCdna.Translate(snp.Chrom == "MT")
		for i := 0; i < cdna.Len(); i++ {
			if j, na1, na2 := i/3, cdna.Base(i), newCdna.Base(i); na1 != na2 && j < protein.Len() {
				aa1, aa2 := protein.Base(j), newProtein.Base(j)
				if aa1 == aa2 {
					anno.Event = "snp_synonymous"
				} else {
					if aa1 == '*' {
						anno.Event = "snp_stoploss"
					} else if aa2 == '*' {
						anno.Event = "snp_stopgain"
					} else {
						if j == 0 {
							anno.Event = "snp_startloss"
						}
						anno.Event = "snp_nonsynonymous"
					}
				}
				anno.NaChange = fmt.Sprintf("c.%d%c>%c", i+1, na1, na2)
				anno.AaChange = fmt.Sprintf("p.%s%d%s", seq.AAMap[aa1], j+1, seq.AAMap[aa2])
				break
			}
		}
	}
	return anno
}

func FindRegion(pos int, regions gene.Regions, strand byte) (index int, exonLen int) {
	for j := 0; j < regions.Len(); j++ {
		index = j
		if strand == '-' {
			index = regions.Len() - j - 1
		}
		region := regions[index]
		if region.Start <= pos && pos <= region.End {
			break
		}
		if region.IsCDS() {
			exonLen += region.End - region.Start + 1
		}
	}
	return index, exonLen
}

func NewAnnotationOfSnp(snp Snv, refgene gene.Refgene) Annotation {
	var anno Annotation
	regionIndex, exonLen := FindRegion(snp.Start, refgene.Regions, refgene.Strand)
	region := refgene.Regions[regionIndex]
	if region.Type == "intron" {
		anno = NewAnnotationOfIntronSnp(snp, refgene, regionIndex, exonLen)
	} else if strings.HasPrefix(region.Type, "UTR") {
		anno = NewAnnotationOfUtrSnp(refgene, regionIndex)
	} else {
		anno = NewAnnotationOfCdsSnp(snp, refgene, regionIndex, exonLen)
	}
	anno.GeneSymbol = refgene.Gene
	anno.GeneEntrezId = refgene.EntrezId
	anno.Transcript = refgene.Transcript
	return anno
}
