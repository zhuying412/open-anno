package snv

import (
	"fmt"
	"grandanno/gene"
	"grandanno/seq"
	"sort"
	"strings"
)

func NewAnnotationOfIntronDel(del Snv, refgene gene.Refgene, regionIndex int, exonLen int) (anno Annotation) {
	region := refgene.Regions[regionIndex]
	prevRegion, _ := refgene.Regions.GetPrev(regionIndex, refgene.Strand)
	nextRegion, _ := refgene.Regions.GetNext(regionIndex, refgene.Strand)
	startDistance1, endDistance1 := del.Start-region.Start+1, del.End-region.Start+1
	startDistance2, endDistance2 := region.End-del.Start+1, region.End-del.End+1
	startDistance, endDistance := startDistance1, endDistance1
	if startDistance1 > endDistance2 {
		startDistance, endDistance = startDistance2, endDistance2
	}
	closestRegion, pos, flag := prevRegion, exonLen, '+'
	if (refgene.Strand == '+') != (startDistance1 <= endDistance2) {
		closestRegion, pos, flag = nextRegion, exonLen+1, '-'
	}
	anno.Region = "intronic"
	if startDistance <= SplicingDistance || endDistance <= SplicingDistance {
		anno.Region = "splicing"
		if !closestRegion.IsCDS() {
			anno.Region = closestRegion.Type + "_splicing"
		}
	}
	if closestRegion.IsCDS() {
		anno.SetExon(closestRegion.ExonOrder)
		anno.NaChange = fmt.Sprintf("c.%d%c%ddel", pos, startDistance, flag)
		if startDistance != endDistance {
			dis1, dis2 := startDistance, endDistance
			if (flag == '+' && dis1 > dis2) || (flag == '-' && dis1 < dis2) {
				dis1, dis2 = endDistance, startDistance
			}
			anno.NaChange = fmt.Sprintf("c.%d%c%d_%d%c%ddel", pos, flag, dis1, flag, dis2)
		}
	}
	return anno
}

func NewAnnotationOfMultiRegionDel(regions gene.Regions) (anno Annotation) {
	regionTypes := make([]string, 0)
	for _, region := range regions {
		if sort.SearchStrings(regionTypes, region.Type) != -1 {
			regionTypes = append(regionTypes, region.Type)
		}
	}
	sort.Strings(regionTypes)
	anno.Region = strings.Join(regionTypes, "_")
	return anno
}

func NewAnnotationOfCdsDel(del Snv, refgene gene.Refgene, lenL int, lenR int, regions gene.Regions) (anno Annotation) {
	lenL, lenR = GetExonLen(lenL, lenR, del.Start, del.End, regions, refgene.Strand)
	cdna, protein := refgene.Cdna, refgene.Protein
	newCdna := cdna.ChangeWithDel(lenL, lenR)
	newProtein := newCdna.Translate(del.Chrom == "MT")
	lenl, lenr := 0, 0
	for i := 0; i < lenL+lenR; i++ {
		if cdna[i] != newCdna[i] {
			break
		}
		lenl++
	}
	for i := lenL + lenR - 1; i >= 0; i-- {
		if cdna[i] != newCdna[i] {
			break
		}
		lenr++
	}
	if lenl+lenr > lenL+lenR {
		lenr = lenL + lenR - lenl
	}
	anno.AaChange = fmt.Sprintf("c.%d_%ddel", lenl+1, lenL+lenR-lenr)
	lenDel := lenL + lenR - lenl - lenr
	lenl, lenr, lenp, lenvp := 0, 0, protein.Len(), newProtein.Len()
	for i := 0; i < lenvp; i++ {
		if protein[i] == newProtein[i] {
			break
		}
		lenl++
	}
	if lenDel%3 == 0 {
		for i := lenvp - 1; i >= 0; i-- {
			if protein[i] != newProtein[i] {
				break
			}
			lenr++
		}
		if lenl+lenr > lenvp {
			lenr = lenvp - lenl
		}
		start, end1, end2 := lenl+1, lenp-lenr, lenvp-lenr
		if newProtein.Find('*') < 0 {
			anno.Event = "del_nonframeshift_stoploss"
		} else if newProtein.Find('*') < lenvp-1 {
			anno.Event = "del_nonframeshift_stopgain"
		} else if protein.Base(0) != protein.Base(0) {
			anno.Event = "del_nonframeshift_startloss"
		} else {
			anno.Event = "del_nonframeshift"
		}
		if start == end2+1 {
			if start == end1 {
				anno.AaChange = fmt.Sprintf("p.%s%ddel", seq.AAMap[protein.Base(start-1)], start)
			} else {
				anno.AaChange = fmt.Sprintf(
					"p.%s%d_%s%ddel",
					seq.AAMap[protein.Base(start-1)],
					start,
					seq.AAMap[protein.Base(end1-1)],
					end1,
				)
			}
		} else {
			anno.AaChange = fmt.Sprintf(
				"p.%s%d_%s%dinsdel%s",
				seq.AAMap[protein.Base(start-1)],
				start,
				seq.AAMap[protein.Base(end1-1)],
				end1,
				newProtein.SubSeq(start-1, end2-start+1).ProteinOne2Tree(),
			)
		}
	} else {
		start := lenl + 1
		if start > newProtein.Len() {
			anno.Event = "del_nonframeshift_stoploss"
			anno.AaChange = fmt.Sprintf(
				"p.%s%d_%s%s%ddel",
				seq.AAMap[protein.Base(start-1)],
				start,
				seq.AAMap[protein[lenp-1]],
				lenp,
			)
		} else {
			if newProtein.Find('*') < 0 {
				anno.Event = "del_frameshift_stoploss"
			} else if newProtein.Find('*') < lenvp-1 {
				anno.Event = "del_frameshift_stopgain"
			} else if protein.Base(0) != protein.Base(0) {
				anno.Event = "del_nonframeshift_startloss"
			} else {
				anno.Event = "del_frameshift"
			}
			anno.AaChange = fmt.Sprintf(
				"p.%s%d%sfs",
				seq.AAMap[protein.Base(start-1)],
				start,
				seq.AAMap[newProtein.Base(start-1)],
			)
		}
	}
	anno.Region = "exonic"
	for _, region := range regions {
		if region.IsCDS() {
			anno.SetExon(region.ExonOrder)
			break
		}
	}
	return anno
}

func FindRegions(start int, end int, regions gene.Regions, strand byte) (indexes []int, lenL int, lenR int) {
	for j := 0; j < regions.Len(); j++ {
		i := j
		if strand == '-' {
			i = regions.Len() - j - 1
		}
		region := regions[i]
		condition1, condition2 := region.Start > end, region.End < start
		if strand == '-' {
			condition1, condition2 = region.End < start, region.Start > end
		}
		if condition1 {
			if region.IsCDS() {
				lenR += region.End - region.Start + 1
			}
		} else if condition2 {
			if region.IsCDS() {
				lenL += region.End - region.Start + 1
			}
		} else {
			indexes = append(indexes, i)
		}
	}
	return indexes, lenL, lenR
}

func GetExonLen(lenL int, lenR int, start int, end int, regions gene.Regions, strand byte) (int, int) {
	for _, region := range regions {
		if region.IsCDS() {
			if region.Start <= start && start <= region.End {
				if strand == '+' {
					lenL += start - region.Start
				} else {
					lenR += start - region.Start
				}
			}
			if region.Start <= end && end <= region.End {
				if strand == '+' {
					lenR += region.End - end
				} else {
					lenL += region.End - end
				}
			}
		}
	}
	return lenL, lenR
}

func NewAnnotationOfDel(del Snv, refgene gene.Refgene) Annotation {
	var anno Annotation
	regionIndexes, lenL, lenR := FindRegions(del.Start, del.End, refgene.Regions, '+')
	regions := make(gene.Regions, 0)
	for _, i := range regionIndexes {
		regions = append(regions, refgene.Regions[i])
	}
	if len(regionIndexes) == 1 {
		regionIndex := regionIndexes[0]
		region := refgene.Regions[regionIndex]
		if region.Type == "intron" {
			anno = NewAnnotationOfIntronDel(del, refgene, regionIndex, lenL)
		} else if strings.HasPrefix(region.Type, "UTR") {
			anno = NewAnnotationOfUtrSnp(refgene, regionIndex)
		} else {
			anno = NewAnnotationOfCdsDel(del, refgene, lenL, lenR, regions)
		}
	} else {
		if regions.HasCds() {
			anno = NewAnnotationOfCdsDel(del, refgene, lenL, lenR, regions)
		} else {
			anno = NewAnnotationOfMultiRegionDel(regions)
		}
	}
	anno.GeneSymbol = refgene.Gene
	anno.GeneEntrezId = refgene.EntrezId
	anno.Transcript = refgene.Transcript
	return anno
}
