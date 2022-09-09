package snv

// import (
// 	"fmt"
// 	"log"
// 	"open-anno/pkg/io"
// 	"open-anno/pkg/schema"
// 	"open-anno/pkg/seq"
// 	"sort"
// 	"strings"
// )

// var AA_SHORT = false

// type TransAnno struct {
// 	Gene       string `json:"gene"`
// 	GeneID     string `json:"gene_id"`
// 	Transcript string `json:"transcript"`
// 	Region     string `json:"region"` // such as: exonic, intronic
// 	NAChange   string `json:"na_change"`
// 	AAChange   string `json:"aa_change"`
// 	Event      string `json:"event"`
// 	Region2    string `json:"region2"` // such as: CDS1, exon1, intron1
// }

// func (this TransAnno) Detail() string {
// 	var detail string
// 	if this.NAChange != "" {
// 		detail = fmt.Sprintf("%s:%s:%s:%s", this.Gene, this.Transcript, this.Region2, this.NAChange)
// 		if this.AAChange != "" {
// 			detail += fmt.Sprintf(":%s", this.AAChange)
// 		}
// 	}
// 	return detail
// }

// func NewTransAnno(trans schema.Transcript, regions ...schema.Region) TransAnno {
// 	anno := TransAnno{
// 		Gene:       trans.Gene,
// 		GeneID:     trans.GeneID,
// 		Transcript: trans.Name,
// 	}
// 	nregions := make(schema.Regions, 0)
// 	for _, region := range regions {
// 		if region.Exists() {
// 			nregions = append(nregions, region)
// 		}
// 	}
// 	if trans.Strand == "+" {
// 		sort.Sort(nregions)
// 	} else {
// 		sort.Sort(sort.Reverse(nregions))
// 	}
// 	if len(nregions) > 0 {
// 		region1, region2 := nregions[0], nregions[len(nregions)-1]
// 		anno.Region2 = "."
// 		if region1.Name() != "" {
// 			anno.Region2 = region1.Name()
// 			if region2.Name() != "" && region1.Name() != region2.Name() {
// 				anno.Region2 = fmt.Sprintf("%s_%s", region1.Name(), region2.Name())
// 			}
// 		} else {
// 			if region2.Name() != "" {
// 				anno.Region2 = region2.Name()
// 			}
// 		}
// 		if region1.Type == schema.RType_CDS || region2.Type == schema.RType_CDS {
// 			anno.Region = "exonic"
// 		} else {
// 			if region1.Type == schema.RType_UTR {
// 				anno.Region = region1.Name()
// 			} else {
// 				if region2.Type == schema.RType_UTR {
// 					anno.Region = region2.Name()
// 				} else {
// 					anno.Region = "intronic"
// 				}
// 			}
// 		}
// 	}
// 	return anno
// }

// type GeneAnno struct {
// 	Gene    string
// 	GeneID  string
// 	Regions []string
// 	Events  []string
// 	Details []string
// }

// func (this *GeneAnno) AddTransAnno(transAnno TransAnno) {
// 	region, event, detail := transAnno.Region, transAnno.Event, transAnno.Detail()
// 	if region == "" {
// 		region = "."
// 	}
// 	if event == "" {
// 		event = "."
// 	}
// 	if detail == "" {
// 		detail = "."
// 	}
// 	this.Regions = append(this.Regions, region)
// 	this.Events = append(this.Events, event)
// 	this.Details = append(this.Details, detail)
// }

// func (this GeneAnno) Region() string {
// 	var regions1, regions2, regions3, regions4 []string
// 	for _, region := range this.Regions {
// 		switch region {
// 		case "exonic", "splicing", "exonic_splicing", "transcript":
// 			regions1 = append(regions1, region)
// 		case "ncRNA":
// 			regions2 = append(regions2, region)
// 		case "UTR3", "UTR5":
// 			regions3 = append(regions3, region)
// 		case "intronic":
// 			regions4 = append(regions4, region)
// 		}
// 	}
// 	if len(regions1) > 0 {
// 		return strings.Join(regions1, ",")
// 	}
// 	if len(regions2) > 0 {
// 		return strings.Join(regions2, ",")
// 	}
// 	if len(regions3) > 0 {
// 		return strings.Join(regions3, ",")
// 	}
// 	if len(regions4) > 0 {
// 		return strings.Join(regions4, ",")
// 	}
// 	return "."
// }

// func (this GeneAnno) Event() string {
// 	if len(this.Events) > 0 {
// 		return strings.Join(this.Events, ",")
// 	}
// 	return "."
// }

// func (this GeneAnno) Detail() string {
// 	if len(this.Details) > 0 {
// 		return strings.Join(this.Details, ",")
// 	}
// 	return "."
// }

// func AnnoSnv(snv schema.Variant, transNames []string, transcripts schema.Transcripts) map[string]GeneAnno {
// 	// var esAnnos, nesAnnos, unkAnnos [TransAnno // exonic_or_splicing, non_exonic_and_non_splicing, ncRNA
// 	var transAnnos []TransAnno
// 	for _, transName := range transNames {
// 		trans := transcripts[transName]
// 		if trans.TxStart <= snv.End && trans.TxEnd >= snv.Start {
// 			var anno TransAnno
// 			if trans.IsUnk() {
// 				anno = NewTransAnno(trans)
// 				anno.Region = "ncRNA"
// 				transAnnos = append(transAnnos, anno)
// 			} else {
// 				if snv.Type() == schema.VType_SNP {
// 					anno = AnnoSnp(snv, trans)
// 				} else if snv.Type() == schema.VType_INS {
// 					anno = AnnoIns(snv, trans)
// 				} else if snv.Type() == schema.VType_DEL {
// 					anno = AnnoDel(snv, trans)
// 				} else {
// 					anno = AnnoSub(snv, trans)
// 				}
// 				transAnnos = append(transAnnos, anno)
// 			}
// 		}
// 	}
// 	geneAnnos := make(map[string]GeneAnno)
// 	for _, transAnno := range transAnnos {
// 		geneAnno, ok := geneAnnos[transAnno.Gene]
// 		if !ok {
// 			geneAnno = GeneAnno{
// 				Gene:   transAnno.Gene,
// 				GeneID: transAnno.GeneID,
// 			}
// 			if geneAnno.GeneID == "" {
// 				geneAnno.GeneID = "."
// 			}
// 		}
// 		geneAnno.AddTransAnno(transAnno)
// 		geneAnnos[transAnno.Gene] = geneAnno
// 	}
// 	return geneAnnos
// }

// func AnnoSnvs(avinput, output, dbname string, geneData GeneData, aashort bool) error {
// 	AA_SHORT = aashort
// 	snvMap, err := io.ReadVariantMap(avinput)
// 	if err != nil {
// 		return err
// 	}
// 	// writer
// 	writer, err := io.NewIoWriter(output)
// 	if err != nil {
// 		return err
// 	}
// 	defer writer.Close()
// 	fmt.Fprintf(writer,
// 		"Chr\tStart\tEnd\tRef\tAlt\t%s.Gene\t%s.GeneID\t%s.Event\t%s.Region\t%s.Detail\n",
// 		dbname, dbname, dbname, dbname, dbname,
// 	)
// 	// anno
// 	for chrom, snvs := range snvMap {
// 		if _, ok := seq.GENOME[chrom]; ok {
// 			log.Printf("Start run annotate %s chr%s ...", dbname, chrom)
// 			if err != nil {
// 				return err
// 			}
// 			transcripts, err := geneData.FilterTranscripts(chrom)
// 			if err != nil {
// 				return err
// 			}
// 			transIndexes := geneData.FilterTransIndexes(chrom)
// 			if err != nil {
// 				return err
// 			}
// 			sort.Sort(snvs)
// 			sort.Sort(transIndexes)
// 			results := make(map[string]map[string]GeneAnno)
// 			for i, j := 0, 0; i < len(snvs) && j < len(transIndexes); {
// 				if snvs[i].End < transIndexes[j].Start {
// 					i++
// 				} else if snvs[i].Start > transIndexes[j].End {
// 					j++
// 				} else {
// 					results[snvs[i].ID()] = AnnoSnv(snvs[i], transIndexes[j].Transcripts, transcripts)
// 					i++
// 				}
// 			}
// 			for _, snv := range snvs {
// 				if geneAnnos, ok := results[snv.ID()]; ok {
// 					for _, geneAnno := range geneAnnos {
// 						fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
// 							snv.Chrom, snv.Start, snv.End, snv.Ref, snv.Alt,
// 							geneAnno.Gene, geneAnno.GeneID, geneAnno.Event(), geneAnno.Region(), geneAnno.Detail(),
// 						)
// 					}
// 				} else {
// 					fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t.\t.\t.\t.\t.\n",
// 						snv.Chrom, snv.Start, snv.End, snv.Ref, snv.Alt,
// 					)
// 				}
// 			}
// 		}
// 	}
// 	return err
// }
