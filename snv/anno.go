package snv

import (
	"bytes"
	"grandanno/core"
	"strconv"
	"strings"
)

type Annotation struct {
	Gene       string
	EntrezId   int
	Transcript string
	Exon       string
	NaChange   string
	AaChange   string
	Region     string
	Function   string
}

type Annotations []Annotation

func (anno *Annotation) SetExon(exonOrder int) {
	var buffer bytes.Buffer
	buffer.WriteString("exon")
	buffer.WriteString(strconv.Itoa(exonOrder))
	anno.Exon = buffer.String()
}

func (annos Annotations) IsSpecial() bool {
	for _, anno := range annos {
		if strings.Contains(anno.Region, "splic") || strings.Contains(anno.Region, "exon") {
			return true
		}
	}
	return false
}

func (annos *Annotations) AnnoIntergeic() {
	*annos = append(*annos, Annotation{Region: "intergenic"})
}

func (annos *Annotations) AnnoStream(snv Snv, refgenes core.Refgenes) {
	for _, refgene := range refgenes {
		for _, region := range refgene.Streams {
			if snv.Variant.Start <= region.End && snv.Variant.End >= region.Start {
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

func (annos *Annotations) AnnoGene(snv Snv, refgenes core.Refgenes, splicingLen int) {
	var cmplAnnos, incmplAnnos, unkAnnos Annotations
	for _, refgene := range refgenes {
		if snv.Variant.End >= refgene.Position.ExonStart && snv.Variant.Start <= refgene.Position.ExonEnd {
			anno := Annotation{
				Gene:       refgene.Gene,
				EntrezId:   refgene.EntrezId,
				Transcript: refgene.Transcript,
			}
			if refgene.Tag == "unk" {
				anno.Region = "unkCDS"
				unkAnnos = append(unkAnnos, anno)
			} else {
				switch snv.GetTypo() {
				case "del":
					anno.AnnoDel(snv, refgene, splicingLen)
				case "ins":
					anno.AnnoIns(snv, refgene, splicingLen)
				case "snp":
					anno.AnnoSnp(snv, refgene, splicingLen)
				default:
					anno.AnnoSnp(snv, refgene, splicingLen)
				}
				if refgene.IsCmpl() {
					cmplAnnos = append(cmplAnnos, anno)
				} else {
					anno.Function = "incmplCDS"
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
