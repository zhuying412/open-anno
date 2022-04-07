package genebased

import (
	"fmt"
	"open-anno/pkg/gene"
	"open-anno/pkg/variant"
	"os"
	"sort"
	"strings"
)

func AnnoCnv(cnv variant.Variant, trans gene.Transcript) CnvGeneBased {
	var cdss, utr3s, utr5s gene.Regions
	var cdsCount int
	regions := trans.Regions
	if trans.Strand == "-" {
		sort.Sort(sort.Reverse(regions))
	}
	for _, region := range regions {
		if region.Type == gene.RType_CDS {
			cdsCount++
		}
		if cnv.Start <= region.End && cnv.End >= region.Start {
			if region.Type == gene.RType_CDS {
				cdss = append(cdss, region)
			}
			if region.Type == gene.RType_UTR {
				if region.Order == 3 {
					utr3s = append(utr3s, region)
				} else {
					utr5s = append(utr5s, region)
				}
			}
		}
	}
	anno := NewCnvGeneBased(trans)
	if len(cdss) > 0 {
		if len(utr5s) > 0 {
			anno.Region = "UTR5_exonic"
			if len(utr3s) > 0 {
				anno.Region = "CDNA"
				if cnv.Start <= trans.TxStart && cnv.End >= trans.TxEnd {
					anno.Region = "transcript"
				}
			}
		} else {
			anno.Region = "exonic"
			if len(utr3s) > 0 {
				anno.Region = "exonic_UTR3"
			}
		}
		if len(cdss) == 1 {
			anno.CDS = fmt.Sprintf("CDS%d/%d", cdss[0].Order, cdsCount)
		} else {
			anno.CDS = fmt.Sprintf("CDS%d_%d/%d", cdss[0].Order, cdss[len(cdss)-1].Order, cdsCount)
		}
	} else {
		if len(utr5s) > 0 {
			anno.Region = "UTR5"
			if len(utr3s) > 0 {
				anno.Region = "ncRNA"
			}
		} else {
			anno.Region = "intronic"
			if len(utr3s) > 0 {
				anno.Region = "UTR3"
			}
		}
	}
	return anno
}

func AnnoCnvs(cnvs variant.Variants, transcripts gene.Transcripts, transIndexes gene.TransIndexes, writer *os.File) {
	sort.Sort(cnvs)
	sort.Sort(transIndexes)
	for _, cnv := range cnvs {
		transNames := make(map[string]bool)
		annos := make([]CnvGeneBased, 0)
		for _, index := range transIndexes {
			if cnv.Start <= index.End && cnv.End >= index.Start {
				for _, transName := range index.Transcripts {
					trans := transcripts[transName]
					if _, ok := transNames[transName]; ok || trans.IsUnk() {
						continue
					}
					transNames[transName] = true
					if trans.IsCmpl() {
						if cnv.Start <= trans.TxEnd && cnv.End >= trans.TxStart {
							anno := AnnoCnv(cnv, trans)
							annos = append(annos, anno)
						}
					}
				}
			}
		}
		annoTexts := make([]string, 0)
		for _, anno := range annos {
			annoTexts = append(annoTexts, fmt.Sprintf("%s:%s:%s:%s:%s:%s:%s",
				anno.Gene, anno.GeneID, anno.Transcript, anno.Strand, anno.Region, anno.CDS, anno.Position,
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
