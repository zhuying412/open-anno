package genebased

import (
	"fmt"
	"open-anno/pkg/gene"
	"open-anno/pkg/variant"
	"os"
	"sort"
)

type SnvGeneBased struct {
	Gene       string `json:"gene"`
	GeneID     string `json:"gene_id"`
	Transcript string `json:"transcript"`
	Region     string `json:"region"`
	NAChange   string `json:"na_change"`
	AAChange   string `json:"aa_change"`
	Event      string `json:"event"`
	Region2    string `json:"region2"`
}

func NewSnvGeneBased(trans gene.Transcript, region gene.Region) SnvGeneBased {
	anno := SnvGeneBased{
		Gene:       trans.Gene,
		GeneID:     trans.GeneID,
		Transcript: trans.Name,
		Event:      ".",
		NAChange:   ".",
		AAChange:   ".",
	}
	anno.Region2 = region.Name()
	if region.Type == gene.RType_UTR {
		anno.Region = region.Name()
	} else if region.Type == gene.RType_INTRON {
		anno.Region = "intronic"
	} else {
		anno.Region = "exonic"
	}
	return anno
}

func fitlerSpecialAnno(annos []SnvGeneBased) []SnvGeneBased {
	filtered := make([]SnvGeneBased, 0)
	for _, anno := range annos {
		if anno.Event != "." {
			filtered = append(filtered, anno)
		}
	}
	return filtered
}

func AnnoSnvs(snvs variant.Variants, transcripts gene.Transcripts, transIndexes gene.TransIndexes, aashort bool, writer *os.File) {
	sort.Sort(snvs)
	sort.Sort(transIndexes)
	for i, j := 0, 0; i < len(snvs) && j < len(transIndexes); {
		if snvs[i].End < transIndexes[j].Start {
			i++
		} else if snvs[i].Start > transIndexes[j].End {
			j++
		} else {
			if snvs[i].Type() != variant.VType_INS {
				i++
				continue
			}
			var cmplAnnoS, cmplAnnos, incmplAnnos, unkAnnos []SnvGeneBased
			for _, transName := range transIndexes[j].Transcripts {
				trans := transcripts[transName]
				if trans.TxStart <= snvs[i].End && trans.TxEnd >= snvs[i].Start {
					var anno SnvGeneBased
					if trans.IsUnk() {
						anno.Region = "ncRNA"
						unkAnnos = append(unkAnnos, anno)
					} else {
						if snvs[i].Type() == variant.VType_SNP {
							anno = AnnoSnp(snvs[i], trans, aashort)
						} else if snvs[i].Type() == variant.VType_INS {
							anno = AnnoIns(snvs[i], transcripts[transName], aashort)
						} else {
							anno = AnnoDel(snvs[i], transcripts[transName], aashort)
						}
						if trans.IsCmpl() {
							if anno.Event == "." {
								cmplAnnos = append(cmplAnnos, anno)
							} else {
								cmplAnnoS = append(cmplAnnoS, anno)
							}
						} else {
							incmplAnnos = append(incmplAnnos, anno)
						}
					}
				}
			}
			annos := cmplAnnoS
			if len(annos) == 0 {
				annos = append(annos, cmplAnnos...)
			}
			if len(annos) == 0 {
				annos = append(annos, incmplAnnos...)
			}
			if len(annos) == 0 {
				annos = append(annos, unkAnnos...)
			}
			if len(annos) == 0 {
				annos = []SnvGeneBased{{
					Gene: ".", GeneID: ".",
					Transcript: ".", Region: ".",
					NAChange: ".", AAChange: ".",
					Event: ".", Region2: ".",
				}}
			}
			for _, anno := range annos {
				fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
					snvs[i].Chrom, snvs[i].Start, snvs[i].End, snvs[i].Ref, snvs[i].Alt,
					anno.Gene, anno.Transcript, anno.Event, anno.Region, anno.NAChange, anno.AAChange,
				)
			}
			i++
		}
	}
}

type CnvGeneBased struct {
	Gene       string `json:"gene"`
	GeneID     string `json:"gene_id"`
	Transcript string `json:"transcript"`
	Region     string `json:"region"`
	Exon       string `json:"exon"`
}

func NewCnvGeneBased(trans gene.Transcript, region gene.Region) CnvGeneBased {
	anno := CnvGeneBased{
		Gene:       trans.Gene,
		GeneID:     trans.GeneID,
		Transcript: trans.Name,
		Region:     region.Name(),
		Exon:       region.Exon,
	}
	return anno
}
