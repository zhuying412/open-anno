package cnv

import (
	"grandanno/gene"
)

type Annotation struct {
	GeneSymbol   string `json:"gene_symbol"`
	GeneEntrezId string `json:"gene_entrez_id"`
	Transcript   string `json:"transcript"`
	Region       string `json:"region"`
	StartExon    int    `json:"StartExon"`
	EndExon      int    `json:"end_exon"`
}

func (a *Annotation) SetExon(exon int) {
	if exon < a.StartExon {
		a.StartExon = exon
	}
	if exon > a.EndExon {
		a.EndExon = exon
	}
}

type Annotations []Annotation

func (a Annotations) HasCoding() bool {
	for _, anno := range a {
		if anno.StartExon > 0 && anno.EndExon > 0 {
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
			if refgene.IsUnk() {
				anno.Region = "ncRNA"
				unkAnnos = append(unkAnnos, anno)
			} else {
				for _, region := range refgene.Regions {
					if cnv.Start <= region.End && cnv.End >= region.Start {
						if region.IsCDS() {
							anno.SetExon(region.ExonOrder)
						}
					}
				}
				if refgene.IsCmpl() {
					anno.Region = "Gene"
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
