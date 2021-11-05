package snv

import (
	"fmt"
	"grandanno/gene"
	"grandanno/input"
	"grandanno/seq"
	"strings"
)

func NewAnnotationOfCnsIns(ins input.Snv, refgene gene.Refgene, regionIndex int, exonLen int) (anno GeneAnno) {
	region := refgene.Regions[regionIndex]
	var pos int
	anno.SetExon(region.ExonOrder)
	anno.Region = "exonic"
	if refgene.Strand == '+' {
		pos = exonLen + ins.Start - region.Start + 1
	} else {
		pos = exonLen + region.End - ins.Start
	}
	if refgene.IsCmpl() {
		cdna, protein := refgene.Cdna, refgene.Protein
		newCdna := cdna.ChangeWithIns(pos, ins.Alt)
		newProtein := newCdna.Translate(ins.Chrom == "MT")

		for i := 0; i < cdna.Len(); i++ {
			if cdna.Base(i) != newCdna.Base(i) {
				anno.AaChange = fmt.Sprintf("c.%dins%s", i, newCdna.SubSeq(i, ins.Alt.Len()))
				break
			}
		}
		lenL, lenR := 0, 0
		for i := 0; i < protein.Len(); i++ {
			if protein.Base(i) != newProtein.Base(i) {
				break
			}
			lenL++
		}
		if lenL < protein.Len() {
			if ins.Alt.Len()%3 == 0 {
				if newProtein.Find('*') < 0 {
					anno.Event = "ins_nonframeshift_stoploss"
				} else if newProtein.Find('*') < protein.Len()-1 {
					anno.Event = "ins_nonframeshift_stopgain"
				} else if protein.Base(0) != newProtein.Base(0) {
					anno.Event = "ins_nonframeshift_startloss"
				} else {
					anno.Event = "ins_nonframeshift"
				}
				for i := protein.Len() - 1; i >= 0; i-- {
					if protein.Base(i) != newProtein.Base(i) {
						break
					}
					lenR++
				}
				if lenL+lenR > protein.Len() {
					lenR = protein.Len() - lenL
				}
				start, end, varEnd := lenL, protein.Len()-lenR, newProtein.Len()-lenR
				altAa := newProtein.SubSeq(start, varEnd-start)
				if start == end {
					anno.AaChange = fmt.Sprintf(
						"p.%s%d_%s%dins%s",
						seq.AAMap[protein.Base(start-1)],
						start,
						seq.AAMap[protein.Base(start)],
						start+1,
						altAa.ProteinOne2Tree(),
					)
				} else {
					refAa := protein.SubSeq(start, end-start)
					anno.AaChange = fmt.Sprintf("p.%s%ddelins%s", refAa.ProteinOne2Tree(), start, altAa.ProteinOne2Tree())
				}
			} else {
				if newProtein.Find('*') < 0 {
					anno.Event = "ins_frameshift_stoploss"
				} else if newProtein.Find('*') < protein.Len()-1 {
					anno.Event = "ins_frameshift_stopgain"
				} else if protein.Base(0) != newProtein.Base(0) {
					anno.Event = "ins_frameshift_startloss"
				} else {
					anno.Event = "ins_frameshift"
				}
				start := lenL
				anno.AaChange = fmt.Sprintf(
					"p.%s%d%sfs",
					seq.AAMap[protein.Base(start)],
					start+1,
					seq.AAMap[newProtein.Base(start)],
				)
			}
		}
	}
	return anno
}

func NewAnnotationOfIns(ins input.Snv, refgene gene.Refgene) GeneAnno {
	var anno GeneAnno
	regionIndex, exonLen := FindRegion(ins.Start, refgene.Regions, refgene.Strand)
	region := refgene.Regions[regionIndex]
	if region.Type == "intron" {
		anno = NewAnnotationOfIntronSnp(ins, refgene, regionIndex, exonLen)
	} else if strings.HasPrefix(region.Type, "UTR") {
		anno = NewAnnotationOfUtrSnp(refgene, regionIndex)
	} else {
		anno = NewAnnotationOfCnsIns(ins, refgene, regionIndex, exonLen)
	}
	anno.GeneSymbol = refgene.Gene
	anno.GeneEntrezId = refgene.EntrezId
	anno.Transcript = refgene.Transcript
	return anno
}
