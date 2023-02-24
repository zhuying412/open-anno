package gene

import (
	"fmt"
	"log"
	"open-anno/anno/variant"
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

func AnnoSnv(snv variant.AnnoVariant, transNames []string, transcripts pkg.Transcripts) []TransAnno {
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
				if snv.Type() == variant.VType_SNP {
					transAnno = AnnoSnp(snv, trans)
				} else if snv.Type() == variant.VType_INS {
					transAnno = AnnoIns(snv, trans)
				} else if snv.Type() == variant.VType_DEL {
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

type GeneAnno struct {
	Gene    string
	GeneID  string
	Regions []string
	Events  []string
	Details []string
}

func (this *GeneAnno) AddEvent(event string) {
	if event != "" && event != "." && pkg.FindArr(this.Events, event) < 0 {
		this.Events = append(this.Events, event)
	}
}
func (this *GeneAnno) AddRegion(region string) {
	if region != "" && region != "." && pkg.FindArr(this.Regions, region) < 0 {
		this.Regions = append(this.Regions, region)
	}
}
func (this *GeneAnno) AddDetail(detail string) {
	if detail != "" && detail != "." && pkg.FindArr(this.Details, detail) < 0 {
		this.Details = append(this.Details, detail)
	}
}

func (this GeneAnno) Region() string {
	var regions1, regions2 []string
	for _, region := range this.Regions {
		switch region {
		case "exonic", "splicing", "exonic_splicing", "transcript":
			regions1 = append(regions1, region)
		case "ncRNA", "UTR3", "UTR5", "intronic":
			regions2 = append(regions2, region)
		}
	}
	if len(regions1) > 0 {
		return strings.Join(regions1, ",")
	}
	if len(regions2) > 0 {
		return strings.Join(regions2, ",")
	}
	return "."
}

func (this GeneAnno) Event() string {
	if len(this.Events) > 0 {
		return strings.Join(this.Events, ",")
	}
	return "."
}

func (this GeneAnno) Detail() string {
	if len(this.Details) > 0 {
		return strings.Join(this.Details, ",")
	}
	return "."
}

func AnnoSnvs(
	snvMap *map[string]variant.SNVs,
	gpes pkg.GenePreds,
	allTransIndexes pkg.TransIndexes,
	genome *faidx.Faidx,
	geneSymbolToID map[string]map[string]string,
	aashort bool) error {
	AA_SHORT = aashort
	// 开始注释
	for chrom, snvs := range *snvMap {
		// if chrom != "chr10" {
		// 	continue
		// }
		log.Printf("Start run annotate GenePred %s ...", chrom)
		transcripts, err := pkg.NewTranscriptsWithSeq(gpes, chrom, geneSymbolToID, genome)
		if err != nil {
			return err
		}
		transIndexes := allTransIndexes.FilterChrom(chrom)
		sort.Sort(snvs)
		sort.Sort(transIndexes)
		for i, j := 0, 0; i < len(snvs) && j < len(transIndexes); {
			if snvs[i].AnnoVariant.End < transIndexes[j].Start {
				i++
			} else if snvs[i].AnnoVariant.Start > transIndexes[j].End {
				j++
			} else {
				transAnnos := AnnoSnv(snvs[i].AnnoVariant, transIndexes[j].Transcripts, transcripts)
				geneAnnos := make(map[string]GeneAnno)
				for _, transAnno := range transAnnos {
					geneAnno, ok := geneAnnos[transAnno.Gene]
					if !ok {
						geneAnno = GeneAnno{Gene: transAnno.Gene, GeneID: transAnno.GeneID}
					}
					geneAnno.AddEvent(transAnno.Event)
					geneAnno.AddRegion(transAnno.Region)
					geneAnno.AddDetail(transAnno.Detail())
					geneAnnos[transAnno.Gene] = geneAnno
				}
				var genes, geneIds, events, regions, details []string
				for _, geneAnno := range geneAnnos {
					genes = append(genes, geneAnno.Gene)
					geneIds = append(geneIds, geneAnno.GeneID)
					events = append(events, geneAnno.Event())
					regions = append(regions, geneAnno.Region())
					details = append(details, geneAnno.Detail())
				}
				snvs[i].Info().Set("GENE", strings.Join(genes, "|"))
				snvs[i].Info().Set("GENE_ID", strings.Join(geneIds, "|"))
				snvs[i].Info().Set("EVENT", strings.Join(events, "|"))
				snvs[i].Info().Set("DETAIL", strings.Join(details, "|"))
				i++
			}
		}
		(*snvMap)[chrom] = snvs
	}
	return nil
}
