package cnv

import (
	"fmt"
	"grandanno/core"
	"strings"
)

type Annotation struct {
	Gene       string
	EntrezId   int
	Transcript string
	Region     string
	Function   string
	Exons      []int
}

func (anno *Annotation) AddExon(exon int) {
	anno.Exons = append(anno.Exons, exon)
}

func (anno *Annotation) GetCds() string {
	switch len(anno.Exons) {
	case 0:
		return "."
	case 1:
		return fmt.Sprintf("exon.%d", anno.Exons[0])
	default:
		min, max := anno.Exons[0], anno.Exons[0]
		for _, exon := range anno.Exons {
			if min > exon {
				min = exon
			}
			if max < exon {
				max = exon
			}
		}
		return fmt.Sprintf("exon.%d_%d", min, max)
	}
}

type Annotations []Annotation

func (annos Annotations) IsSpecial() bool {
	for _, anno := range annos {
		if strings.Contains(anno.Region, "exon") {
			return true
		}
	}
	return false
}

func (annos *Annotations) AnnoIntergeic() {
	*annos = append(*annos, Annotation{Region: "intergenic"})
}

func (annos *Annotations) AnnoStream(cnv Cnv, refgenes core.Refgenes) {
	for _, refgene := range refgenes {
		for _, region := range refgene.Streams {
			if cnv.GetVariant().Start <= region.End && cnv.GetVariant().End >= region.Start {
				*annos = append(*annos, Annotation{
					Gene:       refgene.Gene,
					EntrezId:   refgene.EntrezId,
					Transcript: refgene.Transcript,
					Region:     region.Typo,
				})
				break
			}
		}

	}
}

func (annos *Annotations) AnnoGene(cnv Cnv, refgenes core.Refgenes) {
	var cmplAnnos, incmplAnnos, unkAnnos Annotations
	variant := cnv.GetVariant()
	for _, refgene := range refgenes {
		if variant.End >= refgene.Position.ExonStart && variant.Start <= refgene.Position.ExonEnd {
			anno := Annotation{
				Gene:       refgene.Gene,
				EntrezId:   refgene.EntrezId,
				Transcript: refgene.Transcript,
			}
			if cnv.GetTypo() == "DEL" {
				anno.Function = "Deletion"
			} else {
				anno.Function = "Duplication"
			}
			if refgene.Tag == "unk" {
				anno.Region = "unkCDS"
				unkAnnos = append(unkAnnos, anno)
			} else {
				for _, region := range refgene.Regions {
					if variant.Start <= region.End && variant.End >= region.Start {
						anno.Region = region.Typo
						if region.Typo == "cds" && refgene.Tag == "cmpl" {
							anno.AddExon(region.ExonOrder)
						}
					}
				}
				if refgene.Tag == "cmpl" {
					cmplAnnos = append(cmplAnnos, anno)
				} else {
					anno.Region = "incmplCDS"
					incmplAnnos = append(incmplAnnos, anno)
				}
			}
		}
	}
	*annos = append(*annos, cmplAnnos...)
	if !(*annos).IsSpecial() {
		*annos = append(*annos, incmplAnnos...)
	}
	if len(*annos) == 0 {
		*annos = append(*annos, unkAnnos...)
	}
}
