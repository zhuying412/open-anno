package snv

import (
	"bytes"
	"grandanno/gene"
	"log"
	"strconv"
	"strings"
)

type Annotation struct {
	GeneSymbol   string `json:"gene_symbol"`
	GeneEntrezId string `json:"gene_entrez_id"`
	Transcript   string `json:"transcript"`
	Exon         string `json:"exon"`
	NaChange     string `json:"na_change"`
	AaChange     string `json:"aa_change"`
	Region       string `json:"region"`
	Event        string `json:"event"`
}

func (a *Annotation) SetExon(exonOrder int) {
	var buffer bytes.Buffer
	buffer.WriteString("exon")
	buffer.WriteString(strconv.Itoa(exonOrder))
	a.Exon = buffer.String()
}

func (a Annotation) InCodingOrSplicing() bool {
	return strings.Contains(a.Region, "splic") || strings.Contains(a.Region, "exon")
}

type Annotations []Annotation

func (a Annotations) HasCodingorSplicing() bool {
	for _, anno := range a {
		if anno.InCodingOrSplicing() {
			return true
		}
	}
	return false
}

func NewAnnotationsInIntergeic() Annotations {
	return Annotations{Annotation{Region: "intergenic"}}
}

func NewAnnotationsInUpDownStream(snv Snv, refgenes gene.Refgenes) Annotations {
	annos := make(Annotations, 0)
	for _, refgene := range refgenes {
		for _, region := range refgene.Streams {
			if snv.Start <= region.End && snv.End >= region.Start {
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

func NewAnnotationsInGene(snv Snv, refgenes gene.Refgenes) Annotations {
	var annos, cmplAnnos, incmplAnnos, unkAnnos Annotations
	for _, refgene := range refgenes {
		if snv.End >= refgene.Position.ExonStart && snv.Start <= refgene.Position.ExonEnd {
			var anno Annotation
			if refgene.IsUnk() {
				unkAnnos = append(unkAnnos, Annotation{
					GeneSymbol:   refgene.Gene,
					GeneEntrezId: refgene.EntrezId,
					Transcript:   refgene.Transcript,
					Region:       "unkCDS",
				})
			} else {
				switch snv.Type() {
				case "del":
					anno = NewAnnotationOfDel(snv, refgene)
				case "ins":
					anno = NewAnnotationOfIns(snv, refgene)
				case "snp":
					anno = NewAnnotationOfSnp(snv, refgene)
				default:
					log.Panicf("Unknown SNV type:%s", snv.Type())
				}
				if refgene.IsCmpl() {
					cmplAnnos = append(cmplAnnos, anno)
				} else {
					anno.Event = "incmplCDS"
					incmplAnnos = append(incmplAnnos, anno)
				}
			}
		}
	}
	annos = append(annos, cmplAnnos...)
	if !(annos).HasCodingorSplicing() {
		annos = append(annos, incmplAnnos...)
	}
	if len(annos) == 0 {
		annos = append(annos, unkAnnos...)
	}
	return annos
}
