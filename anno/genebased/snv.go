package genebased

import (
	"fmt"
	"open-anno/pkg/gene"
	"open-anno/pkg/variant"
	"os"
	"sort"
)

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
	fmt.Fprint(writer, "Chr\tStart\tEnd\tRef\tAlt\tGene\tGeneID\tEvent\tRegion\tDetail\n")
	for i, j := 0, 0; i < len(snvs) && j < len(transIndexes); {
		if snvs[i].End < transIndexes[j].Start {
			i++
		} else if snvs[i].Start > transIndexes[j].End {
			j++
		} else {
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
				detail := fmt.Sprintf("%s:%s:%s:%s:%s", anno.Gene, anno.Transcript, anno.Region2, anno.NAChange, anno.AAChange)
				fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
					snvs[i].Chrom, snvs[i].Start, snvs[i].End, snvs[i].Ref, snvs[i].Alt,
					anno.Gene, anno.GeneID, anno.Event, anno.Region, detail,
				)
			}
			i++
		}
	}
}
