package gene

import (
	"fmt"
	"log"
	"open-anno/anno"
	"open-anno/pkg"
	"sort"

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

func AnnoSnv(snv anno.Variant, transNames []string, transcripts pkg.Transcripts) []TransAnno {
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
	snvMap map[string]anno.Variants,
	annoOutputFile string,
	dbName string,
	gpes pkg.GenePreds,
	allTransIndexes pkg.TransIndexes,
	genome *faidx.Faidx,
	geneSymbolToID map[string]map[string]string,
	aashort bool) error {
	AA_SHORT = aashort
	// 打开输出文件
	writer, err := pkg.NewIOWriter(annoOutputFile)
	if err != nil {
		return err
	}
	defer writer.Close()
	// 写入Header
	fmt.Fprintf(writer,
		"Chr\tStart\tEnd\tRef\tAlt\t%s.Gene\t%s.GeneID\t%s.Event\t%s.Region\t%s.Detail\n",
		dbName, dbName, dbName, dbName, dbName,
	)
	// 开始注释
	for chrom, snvs := range snvMap {
		log.Printf("Start run annotate %s %s ...", dbName, chrom)
		if err != nil {
			return err
		}
		transcripts, err := pkg.NewTranscriptsWithSeq(gpes, chrom, geneSymbolToID, genome)
		if err != nil {
			return err
		}
		transIndexes := allTransIndexes.FilterChrom(chrom)
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
	return err
}
