package gene

import (
	"fmt"
	"open-anno/pkg"
	"sort"
	"strings"

	"github.com/brentp/bix"
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

func AnnoCnv(cnv *pkg.CNV, tbx *bix.Bix) (map[string]any, error) {
	annoVar := cnv.AnnoVariant()
	transAnnos := make([]CnvTransAnno, 0)
	query, err := tbx.Query(cnv)
	if err != nil {
		return map[string]any{}, err
	}
	for v, e := query.Next(); e == nil; v, e = query.Next() {
		trans, err := pkg.NewTranscript(fmt.Sprintf("%s", v))
		if err != nil {
			return map[string]any{}, err
		}
		if !trans.IsUnk() {
			trans.SetGeneID()
			err := trans.SetRegions()
			if err != nil {
				return map[string]any{}, err
			}
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
				if annoVar.Start <= region.End && annoVar.End >= region.Start {
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
						if annoVar.Start <= trans.TxStart && annoVar.End >= trans.TxEnd {
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
		}
	}
	query.Close()
	annoTexts := make([]string, 0)
	for _, transAnno := range transAnnos {
		annoTexts = append(annoTexts, fmt.Sprintf("%s:%s:%s:%s:%s:%s:%s",
			transAnno.Gene, transAnno.GeneID, transAnno.Transcript, transAnno.Strand, transAnno.Region, transAnno.CDS, transAnno.Position,
		))
	}
	return map[string]any{"DETAIl": strings.Join(annoTexts, ",")}, nil
}

// func AnnoCnvs(vcfFile string, gpeFile string, goroutines int) (anno.AnnoResult, error) {
// 	annoInfos := make(anno.AnnoInfos)
// 	// 打开句柄
// 	reader, err := pkg.NewIOReader(vcfFile)
// 	if err != nil {
// 		return anno.AnnoResult{}, err
// 	}
// 	defer reader.Close()
// 	vcfReader, err := vcfgo.NewReader(reader, false)
// 	if err != nil {
// 		return anno.AnnoResult{}, err
// 	}
// 	defer vcfReader.Close()
// 	gpeTbx, err := bix.New(gpeFile)
// 	if err != nil {
// 		return anno.AnnoResult{}, err
// 	}
// 	defer gpeTbx.Close()
// 	annoInfosChan := make(chan anno.AnnoInfos, goroutines)
// 	errChan := make(chan error, goroutines)
// 	variants := make([]*pkg.CNV, 0)
// 	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
// 		if len(variant.Chrom()) <= 5 {
// 			variants = append(variants, &pkg.CNV{Variant: *variant})
// 		}
// 	}
// 	multiVariants := pkg.SplitArr(variants, goroutines)
// 	for _, variants := range multiVariants {
// 		go annoCnvs(variants, gpeTbx, annoInfosChan, errChan)
// 	}
// 	for i := 0; i < len(multiVariants); i++ {
// 		err := <-errChan
// 		if err != nil {
// 			return anno.AnnoResult{}, err
// 		}
// 		for pk, annoInfo := range <-annoInfosChan {
// 			annoInfos[pk] = annoInfo
// 		}
// 	}
// 	close(annoInfosChan)
// 	close(errChan)
// 	vcfHeaderInfos := map[string]*vcfgo.Info{
// 		"DETAIL": {
// 			Id:          "DETAIL",
// 			Description: "Gene Detail",
// 			Number:      ".",
// 			Type:        "String",
// 		},
// 	}
// 	return anno.AnnoResult{AnnoInfos: annoInfos, VcfHeaderInfo: vcfHeaderInfos}, err
// }
