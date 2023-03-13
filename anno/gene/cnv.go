package gene

import (
	"fmt"
	"log"
	"open-anno/anno"
	"open-anno/pkg"
	"sort"
	"strings"
)

type CnvTransAnno struct {
	Gene       string `json:"gene"`
	GeneID     string `json:"gene_id"`
	Transcript string `json:"transcript"`
	CDS        string `json:"cds"`
	Region     string `json:"region"`
	Strand     string `json:"strand"`
	Position   string `json:"position"`
}

func NewCnvTransAnno(trans pkg.Transcript) CnvTransAnno {
	transAnno := CnvTransAnno{
		Gene:       trans.Gene,
		GeneID:     trans.GeneID,
		Transcript: trans.Name,
		Strand:     trans.Strand,
		CDS:        ".",
		Region:     ".",
		Position:   fmt.Sprintf("%d-%d", trans.TxStart, trans.TxEnd),
	}
	if transAnno.GeneID == "" {
		transAnno.GeneID = "."
	}
	return transAnno
}

func AnnoCnv(cnv anno.AnnoVariant, trans pkg.Transcript) CnvTransAnno {
	var cdss, utr3s, utr5s pkg.Regions
	var cdsCount int
	regions := trans.Regions
	if trans.Strand == "-" {
		sort.Sort(sort.Reverse(regions))
	}
	for _, region := range regions {
		if region.Type == pkg.RType_CDS {
			cdsCount++
		}
		if cnv.Start <= region.End && cnv.End >= region.Start {
			if region.Type == pkg.RType_CDS {
				cdss = append(cdss, region)
			}
			if region.Type == pkg.RType_UTR {
				if region.Order == 3 {
					utr3s = append(utr3s, region)
				} else {
					utr5s = append(utr5s, region)
				}
			}
		}
	}
	transAnno := NewCnvTransAnno(trans)
	if len(cdss) > 0 {
		if len(utr5s) > 0 {
			transAnno.Region = "UTR5_CDS"
			if len(utr3s) > 0 {
				transAnno.Region = "CDNA"
				if cnv.Start <= trans.TxStart && cnv.End >= trans.TxEnd {
					transAnno.Region = "transcript"
				}
			}
		} else {
			transAnno.Region = "CDS"
			if len(utr3s) > 0 {
				transAnno.Region = "CDS_UTR3"
			}
		}
		if len(cdss) == 1 {
			transAnno.CDS = fmt.Sprintf("CDS%d/%d", cdss[0].Order, cdsCount)
		} else {
			transAnno.CDS = fmt.Sprintf("CDS%d_%d/%d", cdss[0].Order, cdss[len(cdss)-1].Order, cdsCount)
		}
	} else {
		if len(utr5s) > 0 {
			transAnno.Region = "UTR5"
			if len(utr3s) > 0 {
				transAnno.Region = "ncRNA"
			}
		} else {
			transAnno.Region = "intronic"
			if len(utr3s) > 0 {
				transAnno.Region = "UTR3"
			}
		}
	}
	return transAnno
}

func AnnoCnvs(
	variants anno.Variants,
	gpes pkg.GenePreds,
	allTransIndexes pkg.TransIndexes,
	geneSymbolToID map[string]map[string]string) (map[string]map[string]any, error) {
	annoInfos := make(map[string]map[string]any)
	for chrom, cnvs := range variants.AggregateByChrom() {
		log.Printf("Filter GeneBased DB by %s ...", chrom)
		transcripts, err := pkg.NewTranscripts(gpes, chrom, geneSymbolToID)
		if err != nil {
			return annoInfos, err
		}
		transIndexes := allTransIndexes.FilterChrom(chrom)
		sort.Sort(cnvs)
		sort.Sort(transIndexes)
		for _, cnv := range cnvs {
			annoVariant := cnv.AnnoVariant()
			transNames := make(map[string]bool)
			transAnnos := make([]CnvTransAnno, 0)
			for _, index := range transIndexes {
				if annoVariant.Start <= index.End && annoVariant.End >= index.Start {
					for _, transName := range index.Transcripts {
						trans := transcripts[transName]
						if _, ok := transNames[transName]; ok || trans.IsUnk() {
							continue
						}
						transNames[transName] = true
						if !trans.IsUnk() {
							if annoVariant.Start <= trans.TxEnd && annoVariant.End >= trans.TxStart {
								transAnno := AnnoCnv(annoVariant, trans)
								transAnnos = append(transAnnos, transAnno)
							}
						}
					}
				}
			}
			annoTexts := make([]string, 0)
			for _, transAnno := range transAnnos {
				annoTexts = append(annoTexts, fmt.Sprintf("%s:%s:%s:%s:%s:%s:%s",
					transAnno.Gene, transAnno.GeneID, transAnno.Transcript, transAnno.Strand, transAnno.Region, transAnno.CDS, transAnno.Position,
				))
			}
			annoInfos[annoVariant.PK()] = map[string]any{"REGION": strings.Join(annoTexts, ",")}
		}
	}
	return annoInfos, nil
}
