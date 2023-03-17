package gene

import (
	"fmt"
	"open-anno/pkg"
	"sort"
	"strings"

	"github.com/brentp/bix"
	"github.com/brentp/vcfgo"
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

func AnnoCnv(cnv *pkg.Variant, tbx *bix.Bix, pkChan chan string, annoInfoChan chan map[string]any, errChan chan error) {
	annoVar := cnv.AnnoVariant()
	transAnnos := make([]CnvTransAnno, 0)
	query, err := tbx.Query(cnv)
	if err != nil {
		pkChan <- ""
		annoInfoChan <- map[string]any{}
		errChan <- err
		return
	}
	for v, e := query.Next(); e == nil; v, e = query.Next() {
		trans, err := pkg.NewTranscript(fmt.Sprintf("%s", v))
		if err != nil {
			pkChan <- ""
			annoInfoChan <- map[string]any{}
			errChan <- err
			return
		}
		if !trans.IsUnk() && annoVar.Start <= trans.TxEnd && annoVar.End >= trans.TxStart {
			trans.SetGeneID()
			err := trans.SetRegions()
			if err != nil {
				pkChan <- ""
				annoInfoChan <- map[string]any{}
				errChan <- err
				return
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
	annoTexts := make([]string, 0)
	for _, transAnno := range transAnnos {
		annoTexts = append(annoTexts, fmt.Sprintf("%s:%s:%s:%s:%s:%s:%s",
			transAnno.Gene, transAnno.GeneID, transAnno.Transcript, transAnno.Strand, transAnno.Region, transAnno.CDS, transAnno.Position,
		))
	}
	pkChan <- cnv.PK()
	annoInfoChan <- map[string]any{"REGION": strings.Join(annoTexts, ",")}
	errChan <- nil
	return
}

func AnnoCnvs(vcfFile string, gpeFile string, goroutines int) (AnnoInfos, map[string]*vcfgo.Info, error) {
	annoInfos := make(AnnoInfos)
	// 打开句柄
	reader, err := pkg.NewIOReader(vcfFile)
	if err != nil {
		return annoInfos, map[string]*vcfgo.Info{}, err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return annoInfos, map[string]*vcfgo.Info{}, err
	}
	defer vcfReader.Close()
	gpeTbx, err := bix.New(gpeFile)
	if err != nil {
		return annoInfos, map[string]*vcfgo.Info{}, err
	}
	defer gpeTbx.Close()
	pkChan := make(chan string, goroutines)
	annoInfoChan := make(chan map[string]any, goroutines)
	errChan := make(chan error, goroutines)
	size := 0
	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
		go AnnoCnv(&pkg.Variant{Variant: *variant}, gpeTbx, pkChan, annoInfoChan, errChan)
	}
	for i := 0; i < size; i++ {
		err := <-errChan
		if err != nil {
			return AnnoInfos{}, map[string]*vcfgo.Info{}, err
		}
		annoInfos[<-pkChan] = <-annoInfoChan
	}
	close(pkChan)
	close(annoInfoChan)
	close(errChan)
	vcfHeaderInfos := map[string]*vcfgo.Info{
		"REGION": {
			Id:          "REGION",
			Description: "Gene Region",
			Number:      ".",
			Type:        "String",
		},
	}
	return annoInfos, vcfHeaderInfos, nil
}
