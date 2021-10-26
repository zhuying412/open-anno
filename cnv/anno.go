package cnv

import (
	"fmt"
	"grandanno/gene"
	"strings"
)

type Annotation struct {
	GeneSymbol   string `json:"gene_symbol"`
	GeneEntrezId string `json:"gene_entrez_id"`
	Transcript   string `json:"transcript"`
	Region       string `json:"region"`
	Function     string `json:"function"`
	Exons        []int  `json:"exons"`
}

func (a *Annotation) AddExon(exon int) {

	a.Exons = append(a.Exons, exon)
}

func (a Annotation) HasCoding() bool {
	return strings.Contains(a.Region, "exon")
}
func (a *Annotation) GetCds() string {
	switch len(a.Exons) {
	case 0:
		return "."
	case 1:
		return fmt.Sprintf("exon.%d", a.Exons[0])
	default:
		min, max := a.Exons[0], a.Exons[0]
		for _, exon := range a.Exons {
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

func (a Annotations) HasCoding() bool {
	for _, anno := range a {
		if anno.HasCoding() {
			return true
		}
	}
	return false
}

func NewAnnotationsInIntergeic() Annotations {
	return Annotations{Annotation{Region: "intergenic"}}
}

func NewAnnotationsInUpDownStream(cnv Cnv, refgenes gene.Refgenes) Annotations {
	annos := make(Annotations, 0)
	for _, refgene := range refgenes {
		for _, region := range refgene.Streams {
			if cnv.Start <= region.End && cnv.End >= region.Start {
				annos = append(annos, Annotation{
					GeneSymbol:   refgene.Gene,
					GeneEntrezId: refgene.EntrezId,
					Transcript:   refgene.Transcript,
					Region:       region.Type,
				})
				break
			}
		}
	}
	return annos
}

func NewAnnotationsInGene(cnv Cnv, refgenes gene.Refgenes) Annotations {
	var cmplAnnos, incmplAnnos, unkAnnos, annos Annotations
	for _, refgene := range refgenes {
		if cnv.End >= refgene.Position.ExonStart && cnv.Start <= refgene.Position.ExonEnd {
			anno := Annotation{
				GeneSymbol:   refgene.Gene,
				GeneEntrezId: refgene.EntrezId,
				Transcript:   refgene.Transcript,
			}
			if cnv.Ref == "DEL" {
				anno.Function = "Deletion"
			} else {
				anno.Function = "Duplication"
			}
			if refgene.Tag == "unk" {
				anno.Region = "unkCDS"
				unkAnnos = append(unkAnnos, anno)
			} else {
				for _, region := range refgene.Regions {
					if cnv.Start <= region.End && cnv.End >= region.Start {
						anno.Region = region.Type
						if region.Type == "CDS" && refgene.Tag == "cmpl" {
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
	annos = append(annos, cmplAnnos...)
	if !annos.HasCoding() {
		annos = append(annos, incmplAnnos...)
	}
	if len(annos) == 0 {
		annos = append(annos, unkAnnos...)
	}
	return annos
}
