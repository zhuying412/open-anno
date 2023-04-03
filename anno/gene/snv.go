package gene

import (
	"fmt"
	"open-anno/anno"
	"open-anno/pkg"
	"sort"
	"strings"

	"github.com/brentp/bix"
	"github.com/brentp/faidx"
	"github.com/brentp/vcfgo"
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

func annoSnvs(snvs []*pkg.Variant, tbx *bix.Bix, genome *faidx.Faidx, annoInfosChan chan anno.AnnoInfos, errChan chan error) {
	// var esAnnos, nesAnnos, unkAnnos [TransAnno // exonic_or_splicing, non_exonic_and_non_splicing, ncRNA
	annoInfos := make(anno.AnnoInfos)
	for _, snv := range snvs {
		// var esAnnos, nesAnnos, unkAnnos [TransAnno // exonic_or_splicing, non_exonic_and_non_splicing, ncRNA
		annoVar := snv.AnnoVariant()
		query, err := tbx.Query(snv)
		if err != nil {
			annoInfosChan <- anno.AnnoInfos{}
			errChan <- err
			return
		}
		geneAnnos := make(map[string]map[string][]string)
		for v, e := query.Next(); e == nil; v, e = query.Next() {
			trans, err := pkg.NewTranscript(fmt.Sprintf("%s", v))
			if err != nil {
				annoInfosChan <- anno.AnnoInfos{}
				errChan <- err
				return
			}
			if trans.TxStart <= annoVar.End && trans.TxEnd >= annoVar.Start {
				trans.SetGeneID()
				err = trans.SetRegionsWithSeq(genome)
				if err != nil {
					annoInfosChan <- anno.AnnoInfos{}
					errChan <- err
					return
				}
				var transAnno TransAnno
				if trans.IsUnk() {
					transAnno = NewTransAnno(trans)
					transAnno.Region = "ncRNA"
				} else {
					if snv.Type() == pkg.VType_SNP {
						transAnno = AnnoSnp(annoVar, trans)
					} else if snv.Type() == pkg.VType_INS {
						transAnno = AnnoIns(annoVar, trans)
					} else if snv.Type() == pkg.VType_DEL {
						transAnno = AnnoDel(annoVar, trans)
					} else {
						transAnno = AnnoSub(annoVar, trans)
					}
				}
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
		}
		query.Close()
		annoData := make(map[string][]string)
		for _, geneAnno := range geneAnnos {
			for key, val := range geneAnno {
				var value string
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
				if value == "" {
					value = "."
				}
				annoData[key] = append(annoData[key], value)
			}
		}
		annoInfo := make(map[string]any)
		for key, val := range annoData {
			annoInfo[strings.ToUpper(key)] = strings.Join(val, ",")
		}
		annoInfos[snv.PK()] = annoInfo
	}
	annoInfosChan <- annoInfos
	errChan <- nil
	return
}

func AnnoSnvs(vcfFile string, gpeFile string, genome *faidx.Faidx, aashort bool, goroutines int) (anno.AnnoResult, error) {
	AA_SHORT = aashort
	annoInfos := make(anno.AnnoInfos)
	// 打开句柄
	reader, err := pkg.NewIOReader(vcfFile)
	if err != nil {
		return anno.AnnoResult{}, err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return anno.AnnoResult{}, err
	}
	defer vcfReader.Close()
	gpeTbx, err := bix.New(gpeFile)
	if err != nil {
		return anno.AnnoResult{}, err
	}
	defer gpeTbx.Close()
	annoInfosChan := make(chan anno.AnnoInfos, goroutines)
	errChan := make(chan error, goroutines)
	variants := make([]*pkg.Variant, 0)
	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
		if len(variant.Chrom()) <= 5 {
			variants = append(variants, &pkg.Variant{Variant: *variant})
		}
	}
	multiVariants := pkg.SplitArr(variants, goroutines)
	for _, variants := range multiVariants {
		go annoSnvs(variants, gpeTbx, genome, annoInfosChan, errChan)
	}
	for i := 0; i < len(multiVariants); i++ {
		err := <-errChan
		if err != nil {
			return anno.AnnoResult{}, err
		}
		for pk, annoInfo := range <-annoInfosChan {
			annoInfos[pk] = annoInfo
		}
	}
	close(annoInfosChan)
	close(errChan)
	vcfHeaderInfos := map[string]*vcfgo.Info{
		"GENE": {
			Id:          "GENE",
			Description: "Gene Symbol",
			Number:      ".",
			Type:        "String",
		},
		"GENE_ID": {
			Id:          "GENE_ID",
			Description: "Gene Entrez ID",
			Number:      ".",
			Type:        "String",
		},
		"REGION": {
			Id:          "REGION",
			Description: "Region in gene, eg: exonic, intronic, UTR3, UTR5",
			Number:      ".",
			Type:        "String",
		},
		"EVENT": {
			Id:          "EVENT",
			Description: "Variant Event, eg: missense, nonsense, splicing",
			Number:      ".",
			Type:        "String",
		},
		"DETAIL": {
			Id:          "DETAIL",
			Description: "Gene detail, FORMAT=Gene:Transcript:Exon:NA_CHANGE:AA_CHANGE",
			Number:      ".",
			Type:        "String",
		},
	}
	return anno.AnnoResult{AnnoInfos: annoInfos, VcfHeaderInfo: vcfHeaderInfos}, err
}
