package snv

import (
	"bytes"
	"grandanno/db"
	"grandanno/gene"
	"grandanno/input"
	"grandanno/seq"
	"log"
	"strconv"
	"strings"
)

type GeneAnno struct {
	GeneSymbol   string `json:"gene_symbol"`
	GeneEntrezId string `json:"gene_entrez_id"`
	Transcript   string `json:"transcript"`
	Exon         string `json:"exon"`
	NaChange     string `json:"na_change"`
	AaChange     string `json:"aa_change"`
	Region       string `json:"region"`
	Event        string `json:"event"`
}

func (a *GeneAnno) SetExon(exonOrder int) {
	var buffer bytes.Buffer
	buffer.WriteString("exon")
	buffer.WriteString(strconv.Itoa(exonOrder))
	a.Exon = buffer.String()
}

func (a GeneAnno) InCodingOrSplicing() bool {
	return strings.Contains(a.Region, "splic") || strings.Contains(a.Region, "exon")
}

func HasCodingorSplicingRegion(geneAnnos []GeneAnno) bool {
	for _, anno := range geneAnnos {
		if anno.InCodingOrSplicing() {
			return true
		}
	}
	return false
}

func NewGeneAnnosInIntergeic() []GeneAnno {
	return []GeneAnno{{Region: "intergenic"}}
}

func NewGeneAnnosInUpDownStream(snv input.Snv, refgenes []gene.Refgene) []GeneAnno {
	annos := make([]GeneAnno, 0)
	for _, refgene := range refgenes {
		for _, region := range refgene.Streams {
			if snv.Start <= region.End && snv.End >= region.Start {
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

func NewGeneAnnosInGene(snv input.Snv, refgenes []gene.Refgene) []GeneAnno {
	var annos, cmplAnnos, incmplAnnos, unkAnnos []GeneAnno
	for _, refgene := range refgenes {
		if snv.End >= refgene.Position.ExonStart && snv.Start <= refgene.Position.ExonEnd {
			var anno GeneAnno
			if refgene.IsUnk() {
				unkAnnos = append(unkAnnos, GeneAnno{
					GeneSymbol:   refgene.Gene,
					GeneEntrezId: refgene.EntrezId,
					Transcript:   refgene.Transcript,
					Region:       "ncRNA",
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
	if !HasCodingorSplicingRegion(annos) {
		annos = append(annos, incmplAnnos...)
	}
	if len(annos) == 0 {
		annos = append(annos, unkAnnos...)
	}
	return annos
}

type GeneAnnoMap map[string][]GeneAnno

func (g *GeneAnnoMap) RunAnnotate(snvs input.Snvs, refgeneMap gene.RefgeneMap, refIndexes db.RefIndexes) {
	for i, j := 0, 0; i < len(snvs) && j < len(refIndexes); {
		snvStart, snvEnd := snvs[i].Start, snvs[i].End
		refStart, refEnd := refIndexes[j].Start, refIndexes[j].End
		if snvStart > refEnd {
			j++
		} else {
			var annos []GeneAnno
			if snvEnd < refStart {
				annos = NewGeneAnnosInIntergeic()
			} else {
				var refgenes []gene.Refgene
				refgenes = refgeneMap.FindMany(refIndexes[j].Transcripts)
				annos = NewGeneAnnosInGene(snvs[i], refgenes)
				if len(annos) == 0 {
					annos = NewGeneAnnosInUpDownStream(snvs[i], refgenes)
				}
				if len(annos) == 0 {
					annos = NewGeneAnnosInIntergeic()
				}
			}
			(*g)[snvs[i].SN()] = annos
			i++
		}
	}
}

func NewGeneAnnoMap(snvs input.Snvs, refgeneMap gene.RefgeneMap, refIndexes db.RefIndexes) GeneAnnoMap {
	InitSnvParam()
	geneAnnoMap := make(GeneAnnoMap)
	for _, chrom := range db.ChromArray {
		subSnvs := snvs.FilterByChrom(chrom.Name)
		if subSnvs.Len() == 0 {
			continue
		}
		mrnaFastaFile := db.GetMrnaFile(chrom)
		log.Printf("read %s\n", mrnaFastaFile)
		mrna := seq.ReadFastaFile(mrnaFastaFile)
		log.Printf("annotate chr%s\n", chrom.Name)
		subRefgeneMap := refgeneMap.FilterByChrom(chrom.Name)
		subRefIndexes := refIndexes.FilterByChrom(chrom.Name)
		refgeneMap.SetSequence(mrna)
		geneAnnoMap.RunAnnotate(subSnvs, subRefgeneMap, subRefIndexes)
	}
	return geneAnnoMap
}
