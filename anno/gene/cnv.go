package gene

import (
	"fmt"
	"log"
	"open-anno/anno/variant"
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

func AnnoCnv(cnv variant.AnnoVariant, trans pkg.Transcript) CnvTransAnno {
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
	cnvMap map[string]variant.CNVs,
	gpes pkg.GenePreds,
	allTransIndexes pkg.TransIndexes,
	geneSymbolToID map[string]map[string]string,
	annoOutput, dbname string) error {
	// 打开输出文件
	writer, err := pkg.NewIOWriter(annoOutput)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprintf(writer, "Chr\tStart\tEnd\tRef\tAlt\t%s.Region\n", dbname)
	for chrom, cnvs := range cnvMap {
		log.Printf("Filter GeneBased DB by %s ...", chrom)
		transcripts, err := pkg.NewTranscripts(gpes, chrom, geneSymbolToID)
		if err != nil {
			return err
		}
		transIndexes := allTransIndexes.FilterChrom(chrom)
		sort.Sort(cnvs)
		sort.Sort(transIndexes)
		for _, cnv := range cnvs {
			transNames := make(map[string]bool)
			transAnnos := make([]CnvTransAnno, 0)
			for _, index := range transIndexes {
				if cnv.Start <= index.End && cnv.End >= index.Start {
					for _, transName := range index.Transcripts {
						trans := transcripts[transName]
						if _, ok := transNames[transName]; ok || trans.IsUnk() {
							continue
						}
						transNames[transName] = true
						if !trans.IsUnk() {
							if cnv.Start <= trans.TxEnd && cnv.End >= trans.TxStart {
								transAnno := AnnoCnv(cnv, trans)
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
			if len(annoTexts) == 0 {
				annoTexts = []string{"."}
			}
			fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\n",
				cnv.Chrom, cnv.Start, cnv.End, cnv.Ref, cnv.Alt, strings.Join(annoTexts, ","),
			)
		}
	}
	return err
}
