package gene

import (
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"open-anno/pkg/io/refgene"
	"sort"
	"strings"
)

type SnvGeneBasedResult struct {
	Gene    string
	GeneID  string
	Regions []string
	Events  []string
	Details []string
}

func AnnoSnvs(snvs io.Variants, transcripts refgene.Transcripts, transIndexes refgene.TransIndexes, aashort bool, writer io.WriteCloser) {
	sort.Sort(snvs)
	sort.Sort(transIndexes)
	results := make(map[string]map[string]SnvGeneBasedResult)
	for i, j := 0, 0; i < len(snvs) && j < len(transIndexes); {
		if snvs[i].End < transIndexes[j].Start {
			i++
		} else if snvs[i].Start > transIndexes[j].End {
			j++
		} else {
			var esAnnos, nesAnnos, unkAnnos []SnvGeneBased // exonic_splicing, non_exonic_splicing, ncRNA
			for _, transName := range transIndexes[j].Transcripts {
				trans := transcripts[transName]
				if trans.TxStart <= snvs[i].End && trans.TxEnd >= snvs[i].Start {
					var anno SnvGeneBased
					if trans.IsUnk() {
						anno.Region = "ncRNA"
						unkAnnos = append(unkAnnos, anno)
					} else {
						if snvs[i].Type() == io.VType_SNP {
							anno = AnnoSnp(snvs[i], trans, aashort)
						} else if snvs[i].Type() == io.VType_INS {
							anno = AnnoIns(snvs[i], transcripts[transName], aashort)
						} else if snvs[i].Type() == io.VType_DEL {
							anno = AnnoDel(snvs[i], transcripts[transName], aashort)
						} else {
							anno = AnnoSub(snvs[i], transcripts[transName], aashort)
						}
						if anno.Event == "" {
							esAnnos = append(esAnnos, anno)
						} else {
							nesAnnos = append(nesAnnos, anno)
						}
					}
				}
			}
			annos := esAnnos
			if len(annos) == 0 {
				annos = append(annos, unkAnnos...)
			}
			if len(annos) == 0 {
				annos = append(annos, nesAnnos...)
			}
			id := snvs[i].ID()
			for _, anno := range annos {
				detail := "."
				if anno.NAChange != "" {
					detail = fmt.Sprintf("%s:%s:%s:%s", anno.Gene, anno.Transcript, anno.Region2, anno.NAChange)
					if anno.AAChange != "" {
						detail += fmt.Sprintf(":%s", anno.AAChange)
					}
				}
				if _, ok := results[id]; !ok {
					results[id] = make(map[string]SnvGeneBasedResult)
				}
				result, ok := results[id][anno.Gene]
				if !ok {
					result = SnvGeneBasedResult{
						Gene:   anno.Gene,
						GeneID: anno.GeneID,
					}
				}
				if anno.Region != "" && pkg.FindArr(result.Regions, anno.Region) < 0 {
					result.Regions = append(result.Regions, anno.Region)
				}
				if anno.Event != "" && pkg.FindArr(result.Events, anno.Event) < 0 {
					result.Events = append(result.Events, anno.Event)
				}
				if detail != "" && pkg.FindArr(result.Details, detail) < 0 {
					result.Details = append(result.Details, detail)
				}

				results[id][anno.Gene] = result
			}
			i++
		}
	}
	for _, snv := range snvs {
		if subResults, ok := results[snv.ID()]; ok {
			for _, result := range subResults {
				geneId, event, region, detail := ".", ".", ".", "."
				if result.GeneID != "" {
					geneId = result.GeneID
				}
				if len(result.Events) > 0 {
					event = strings.Join(result.Events, ",")
				}
				if len(result.Regions) > 0 {
					region = strings.Join(result.Regions, ",")
				}
				if len(result.Details) > 0 {
					detail = strings.Join(result.Details, ",")
				}
				fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
					snv.Chrom, snv.Start, snv.End, snv.Ref, snv.Alt,
					result.Gene, geneId, event, region, detail,
				)
			}
		} else {
			fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t.\t.\t.\t.\t.\n",
				snv.Chrom, snv.Start, snv.End, snv.Ref, snv.Alt,
			)
		}
	}
}
