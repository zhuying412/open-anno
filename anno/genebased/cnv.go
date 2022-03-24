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
	anno := NewCnvGeneBased(trans)
	if cnv.Start <= trans.TxStart {
		if cnv.End >= trans.TxEnd {
			anno.Region = "transcript"
		} else {
			if trans.Strand == "+" {
				anno.Region = "UTR5"
			} else {
				anno.Region = "UTR3"
			}
		}
	} else {
		if cnv.End >= trans.TxEnd {
			if trans.Strand == "+" {
				anno.Region = "UTR3"
			} else {
				anno.Region = "UTR5"
			}
		}
	}
	var cds1, cds2 gene.Region
	var cdsCount int
	regions := trans.Regions
	if trans.Strand == "-" {
		sort.Sort(sort.Reverse(regions))
	}
	for _, region := range regions {
		if region.Type == gene.RType_CDS {
			cdsCount++
			if cnv.Start <= region.End && cnv.End >= region.Start {
				if !cds1.Exists() {
					cds1 = region
				}
				cds2 = region
			}
		}
	}
	if cds1.Exists() {
		if anno.Region != "transcript" {
			anno.Region = "exonic"
		}
		if cds1.Equal(cds2) {
			anno.CDS = fmt.Sprintf("CDS%d/%d", cds1.Order, cdsCount)
		} else {
			anno.CDS = fmt.Sprintf("CDS%d_%d/%d", cds1.Order, cds2.Order, cdsCount)
		}
	}
	return anno
}

func AnnoCnvs(cnvs variant.Variants, transcripts gene.Transcripts, transIndexes gene.TransIndexes, writer *os.File) {
	sort.Sort(cnvs)
	sort.Sort(transIndexes)
	for _, cnv := range cnvs {
		transNames := make([]string, 0)
		annos := make([]CnvGeneBased, 0)
		for _, index := range transIndexes {
			if cnv.Start <= index.End && cnv.End >= index.Start {
				for _, transName := range index.Transcripts {
					if sort.SearchStrings(transNames, transName) < 0 {
						continue
					}
					trans := transcripts[transName]
					if trans.IsCmpl() {
						if cnv.Start <= trans.TxStart && cnv.End >= trans.TxEnd {
							anno := AnnoCnv(cnv, trans)
							annos = append(annos, anno)
						}
					}
				}
			}
		}
		annoTexts := make([]string, 0)
		for _, anno := range annos {
			annoTexts = append(annoTexts,
				fmt.Sprintf("%s:%s:%s:%s:%s:%s", anno.Gene, anno.GeneID, anno.Transcript, anno.Strand, anno.Region, anno.CDS),
			)
		}
		fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\n",
			cnv.Chrom, cnv.Start, cnv.End, cnv.Ref, cnv.Alt, strings.Join(annoTexts, ","),
		)

	}
}
