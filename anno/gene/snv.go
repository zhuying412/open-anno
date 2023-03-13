package gene

import (
	"fmt"
	"log"
	"open-anno/anno"
	"open-anno/pkg"
	"sort"
	"strings"

	"github.com/brentp/faidx"
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

func NewTransAnno(trans pkg.Transcript, regions ...pkg.Region) TransAnno {
	transAnno := TransAnno{
		Gene:       trans.Gene,
		GeneID:     trans.GeneID,
		Transcript: trans.Name,
	}
	nregions := make(pkg.Regions, 0)
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
		transAnno.Region2 = "."
		if region1.Name() != "" {
			transAnno.Region2 = region1.Name()
			if region2.Name() != "" && region1.Name() != region2.Name() {
				transAnno.Region2 = fmt.Sprintf("%s_%s", region1.Name(), region2.Name())
			}
		} else {
			if region2.Name() != "" {
				transAnno.Region2 = region2.Name()
			}
		}
		if region1.Type == pkg.RType_CDS || region2.Type == pkg.RType_CDS {
			transAnno.Region = "exonic"
		} else {
			if region1.Type == pkg.RType_UTR {
				transAnno.Region = region1.Name()
			} else {
				if region2.Type == pkg.RType_UTR {
					transAnno.Region = region2.Name()
				} else {
					transAnno.Region = "intronic"
				}
			}
		}
	}
	return transAnno
}

func AnnoSnv(snv anno.AnnoVariant, transNames []string, transcripts pkg.Transcripts) []TransAnno {
	// var esAnnos, nesAnnos, unkAnnos [TransAnno // exonic_or_splicing, non_exonic_and_non_splicing, ncRNA
	var transAnnos []TransAnno
	for _, transName := range transNames {
		trans := transcripts[transName]
		if trans.TxStart <= snv.End && trans.TxEnd >= snv.Start {
			var transAnno TransAnno
			if trans.IsUnk() {
				transAnno = NewTransAnno(trans)
				transAnno.Region = "ncRNA"
				transAnnos = append(transAnnos, transAnno)
			} else {
				if snv.Type() == anno.VType_SNP {
					transAnno = AnnoSnp(snv, trans)
				} else if snv.Type() == anno.VType_INS {
					transAnno = AnnoIns(snv, trans)
				} else if snv.Type() == anno.VType_DEL {
					transAnno = AnnoDel(snv, trans)
				} else {
					transAnno = AnnoSub(snv, trans)
				}
				transAnnos = append(transAnnos, transAnno)
			}
		}
	}
	return transAnnos
}

func AnnoSnvs(
	variants anno.Variants,
	gpes pkg.GenePreds,
	allTransIndexes pkg.TransIndexes,
	genome *faidx.Faidx,
	geneSymbolToID map[string]map[string]string,
	aashort bool) (map[string]map[string]any, error) {
	AA_SHORT = aashort
	annoInfos := make(map[string]map[string]any)
	// 开始注释
	for chrom, snvs := range variants.AggregateByChrom() {
		// if chrom != "chr10" {
		// 	continue
		// }
		log.Printf("Start run annotate GenePred %s ...", chrom)
		transcripts, err := pkg.NewTranscriptsWithSeq(gpes, chrom, geneSymbolToID, genome)
		if err != nil {
			return annoInfos, err
		}
		transIndexes := allTransIndexes.FilterChrom(chrom)
		sort.Sort(snvs)
		sort.Sort(transIndexes)
		for i, j := 0, 0; i < len(snvs) && j < len(transIndexes); {
			annoVariant := snvs[i].AnnoVariant()
			if annoVariant.End < transIndexes[j].Start {
				i++
			} else if annoVariant.Start > transIndexes[j].End {
				j++
			} else {
				transAnnos := AnnoSnv(annoVariant, transIndexes[j].Transcripts, transcripts)
				geneAnnos := make(map[string]map[string][]string)
				for _, transAnno := range transAnnos {
					geneAnno, ok := geneAnnos[transAnno.Gene]
					if !ok {
						geneAnno = map[string][]string{"gene": {transAnno.Gene}, "gene_id": {transAnno.GeneID}, "region": {}, "event": {}, "detail": {}}
					}
					region, event, detail := transAnno.Region, transAnno.Event, transAnno.Detail()
					if region != "" && region != "." && pkg.FindArr(geneAnno["region"], region) < 0 {
						geneAnno["region"] = append(geneAnno["region"], region)
					}
					if event != "" && event != "." && pkg.FindArr(geneAnno["event"], event) < 0 {
						geneAnno["event"] = append(geneAnno["event"], event)
					}
					if detail != "" && detail != "." && pkg.FindArr(geneAnno["detail"], detail) < 0 {
						geneAnno["detail"] = append(geneAnno["detail"], detail)
					}
					geneAnnos[transAnno.Gene] = geneAnno
				}
				annoData := make(map[string][]string)
				for _, geneAnno := range geneAnnos {
					for key, val := range geneAnno {
						value := "."
						if key == "region" {
							var regions1, regions2 []string
							for _, region := range val {
								switch region {
								case "exonic", "splicing", "exonic_splicing", "transcript":
									regions1 = append(regions1, region)
								case "ncRNA", "UTR3", "UTR5", "intronic":
									regions2 = append(regions2, region)
								}
							}
							if len(regions1) > 0 {
								value = strings.Join(regions1, "|")
							} else {
								if len(regions2) > 0 {
									value = strings.Join(regions2, "|")
								}
							}
						} else {
							value = strings.Join(val, "|")
						}
						annoData[key] = append(annoData[key], value)
					}
				}
				pk := annoVariant.PK()
				for key, val := range annoData {
					_, ok := annoInfos[pk]
					if !ok {
						annoInfos[pk] = make(map[string]any)
					}
					annoInfos[pk][strings.ToUpper(key)] = strings.Join(val, ",")
				}
				i++
			}
		}
	}
	return annoInfos, nil
}
