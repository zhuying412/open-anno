package snv

import (
	"fmt"
	"log"
	"open-anno/pkg/io"
	"open-anno/pkg/schema"
	"open-anno/pkg/seq"
	"sort"
)

var AA_SHORT = false

type TransAnno struct {
	Gene       string `json:"gene"`
	GeneID     string `json:"gene_id"`
	Transcript string `json:"transcript"`
	Region     string `json:"region"` // such as: exonic, intronic
	NAChange   string `json:"na_change"`
	AAChange   string `json:"aa_change"`
	Event      string `json:"event"`
	Region2    string `json:"region2"` // such as: CDS1, exon1, intron1
}

func (this TransAnno) Detail() string {
	var detail string
	if this.NAChange != "" {
		detail = fmt.Sprintf("%s:%s:%s:%s", this.Gene, this.Transcript, this.Region2, this.NAChange)
		if this.AAChange != "" {
			detail += fmt.Sprintf(":%s", this.AAChange)
		}
	}
	return detail
}

func NewTransAnno(trans schema.Transcript, regions ...schema.Region) TransAnno {
	anno := TransAnno{
		Gene:       trans.Gene,
		GeneID:     trans.GeneID,
		Transcript: trans.Name,
	}
	nregions := make(schema.Regions, 0)
	for _, region := range regions {
		if region.Exists() {
			nregions = append(nregions, region)
		}
	}
	if trans.Strand == "+" {
		sort.Sort(nregions)
	} else {
		sort.Sort(sort.Reverse(nregions))
	}
	if len(nregions) > 0 {
		region1, region2 := nregions[0], nregions[len(nregions)-1]
		anno.Region2 = "."
		if region1.Name() != "" {
			anno.Region2 = region1.Name()
			if region2.Name() != "" && region1.Name() != region2.Name() {
				anno.Region2 = fmt.Sprintf("%s_%s", region1.Name(), region2.Name())
			}
		} else {
			if region2.Name() != "" {
				anno.Region2 = region2.Name()
			}
		}
		if region1.Type == schema.RType_CDS || region2.Type == schema.RType_CDS {
			anno.Region = "exonic"
		} else {
			if region1.Type == schema.RType_UTR {
				anno.Region = region1.Name()
			} else {
				if region2.Type == schema.RType_UTR {
					anno.Region = region2.Name()
				} else {
					anno.Region = "intronic"
				}
			}
		}
	}
	return anno
}

func AnnoSnv(snv schema.Variant, transNames []string, transcripts schema.Transcripts) []TransAnno {
	// var esAnnos, nesAnnos, unkAnnos [TransAnno // exonic_or_splicing, non_exonic_and_non_splicing, ncRNA
	var transAnnos []TransAnno
	for _, transName := range transNames {
		trans := transcripts[transName]
		if trans.TxStart <= snv.End && trans.TxEnd >= snv.Start {
			var anno TransAnno
			if trans.IsUnk() {
				anno = NewTransAnno(trans)
				anno.Region = "ncRNA"
				transAnnos = append(transAnnos, anno)
			} else {
				if snv.Type() == schema.VType_SNP {
					anno = AnnoSnp(snv, trans)
				} else if snv.Type() == schema.VType_INS {
					anno = AnnoIns(snv, trans)
				} else if snv.Type() == schema.VType_DEL {
					anno = AnnoDel(snv, trans)
				} else {
					anno = AnnoSub(snv, trans)
				}
				transAnnos = append(transAnnos, anno)
			}
		}
	}
	return transAnnos
}

func AnnoSnvs(avinput, output, dbname string, geneData GeneData, aashort bool) error {
	AA_SHORT = aashort
	snvMap, err := io.ReadVariantMap(avinput)
	if err != nil {
		return err
	}
	// writer
	writer, err := io.NewIoWriter(output)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprintf(writer,
		"Chr\tStart\tEnd\tRef\tAlt\t%s.Gene\t%s.GeneID\t%s.Event\t%s.Region\t%s.Detail\n",
		dbname, dbname, dbname, dbname, dbname,
	)
	// anno
	for chrom, snvs := range snvMap {
		if _, ok := seq.GENOME[chrom]; ok {
			log.Printf("Start run annotate %s chr%s ...", dbname, chrom)
			if err != nil {
				return err
			}
			transcripts, err := geneData.FilterTranscripts(chrom)
			if err != nil {
				return err
			}
			transIndexes := geneData.FilterTransIndexes(chrom)
			if err != nil {
				return err
			}
			sort.Sort(snvs)
			sort.Sort(transIndexes)
			for i, j := 0, 0; i < len(snvs) && j < len(transIndexes); {
				if snvs[i].End < transIndexes[j].Start {
					i++
				} else if snvs[i].Start > transIndexes[j].End {
					j++
				} else {
					chrom, start, end, ref, alt := snvs[i].Chrom, snvs[i].Start, snvs[i].End, snvs[i].Ref, snvs[i].Alt
					transAnnos := AnnoSnv(snvs[i], transIndexes[j].Transcripts, transcripts)
					for _, transAnno := range transAnnos {
						gene, geneId, event, region, detail := transAnno.Gene, transAnno.GeneID, transAnno.Event, transAnno.Region, transAnno.Detail()
						if geneId == "" {
							geneId = "."
						}
						if event == "" {
							event = "."
						}
						if region == "" {
							region = "."
						}
						if detail == "" {
							detail = "."
						}
						fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
							chrom, start, end, ref, alt, gene, geneId, event, region, detail,
						)
					}
					i++
				}
			}
		}
	}
	return err
}
