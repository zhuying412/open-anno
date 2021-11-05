package cnv

import (
	"grandanno/db"
	"grandanno/gene"
	"grandanno/input"
	"log"
)

type GeneAnno struct {
	GeneSymbol   string `json:"gene_symbol"`
	GeneEntrezId string `json:"gene_entrez_id"`
	Transcript   string `json:"transcript"`
	Region       string `json:"region"`
	StartExon    int    `json:"StartExon"`
	EndExon      int    `json:"end_exon"`
}

func (a *GeneAnno) SetExon(exon int) {
	if exon < a.StartExon {
		a.StartExon = exon
	}
	if exon > a.EndExon {
		a.EndExon = exon
	}
}

func HasCodingRegion(geneAnnos []GeneAnno) bool {
	for _, anno := range geneAnnos {
		if anno.StartExon > 0 && anno.EndExon > 0 {
			return true
		}
	}
	return false
}

func NewAnnotationsInIntergeic() []GeneAnno {
	return []GeneAnno{{Region: "intergenic"}}
}

func NewAnnotationsInUpDownStream(cnv input.Cnv, refgenes []gene.Refgene) []GeneAnno {
	annos := make([]GeneAnno, 0)
	for _, refgene := range refgenes {
		for _, region := range refgene.Streams {
			if cnv.Start <= region.End && cnv.End >= region.Start {
				annos = append(annos, GeneAnno{
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

func NewAnnotationsInGene(cnv input.Cnv, refgenes []gene.Refgene) []GeneAnno {
	var cmplAnnos, incmplAnnos, unkAnnos, annos []GeneAnno
	for _, refgene := range refgenes {
		if cnv.End >= refgene.Position.ExonStart && cnv.Start <= refgene.Position.ExonEnd {
			anno := GeneAnno{
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
	if !HasCodingRegion(annos) {
		annos = append(annos, incmplAnnos...)
	}
	if len(annos) == 0 {
		annos = append(annos, unkAnnos...)
	}
	return annos
}

type GeneAnnoMap map[string][]GeneAnno

func (g *GeneAnnoMap) RunAnnotate(cnvs input.Cnvs, refgeneMap gene.RefgeneMap, refIndexes db.RefIndexes) {
	for i := 0; i < len(cnvs); i++ {
		cnvStart, cnvEnd := cnvs[i].Start, cnvs[i].End
		for j := 0; j < len(refIndexes); j++ {
			refStart, refEnd := refIndexes[j].Start, refIndexes[j].End
			if cnvEnd < refStart {
				break
			} else if cnvStart > refEnd {
				continue
			} else {
				var annos []GeneAnno
				if cnvEnd < refStart {
					annos = NewAnnotationsInIntergeic()
				} else {
					refgenes := refgeneMap.FindMany(refIndexes[j].Transcripts)
					annos = NewAnnotationsInGene(cnvs[i], refgenes)
					if len(annos) == 0 {
						annos = NewAnnotationsInUpDownStream(cnvs[i], refgenes)
					}
					if len(annos) == 0 {
						annos = NewAnnotationsInIntergeic()
					}
				}
				(*g)[cnvs[i].SN()] = annos
			}
		}
	}
}

func NewGeneAnnoMap(cnvs input.Cnvs, refgeneMap gene.RefgeneMap, refIndexes db.RefIndexes) GeneAnnoMap {
	geneAnnoMap := make(GeneAnnoMap)
	for _, chrom := range db.ChromArray {
		subCnvs := cnvs.FilterByChrom(chrom.Name)
		if subCnvs.Len() == 0 {
			continue
		}
		log.Printf("annotate chr%s\n", chrom.Name)
		subRefgeneMap := refgeneMap.FilterByChrom(chrom.Name)
		subRefIndexes := refIndexes.FilterByChrom(chrom.Name)
		geneAnnoMap.RunAnnotate(subCnvs, subRefgeneMap, subRefIndexes)
	}
	return geneAnnoMap
}
