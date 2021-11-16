package snp

import (
	"OpenAnno/db/transcript"
	"OpenAnno/variant"
	"bytes"
	"github.com/spf13/viper"
	"strconv"
	"strings"
)

var SplicingDistance int

func Init() {
	SplicingDistance = viper.GetInt("param.splicing_distance")
}

type GeneAnnoItem struct {
	GeneSymbol   string `json:"gene_symbol"`
	GeneEntrezId string `json:"gene_entrez_id"`
	Transcript   string `json:"transcript"`
	Region       string `json:"region"`
	Exon         string `json:"exon"`
	NAChange     string `json:"na_change"`
	AAChange     string `json:"aa_change"`
	Event        string `json:"event"`
}

func (a *GeneAnnoItem) SetExon(exonOrder int) {
	var buffer bytes.Buffer
	buffer.WriteString("exon")
	buffer.WriteString(strconv.Itoa(exonOrder))
	a.Exon = buffer.String()
}

func (a GeneAnnoItem) InCodingOrSplicingRegion() bool {
	return strings.Contains(a.Region, "splic") || strings.Contains(a.Region, "exon")
}

func (a *GeneAnnoItem) SetGene(trans transcript.Transcript) {
	a.GeneSymbol = trans.Gene
	a.GeneEntrezId = trans.EntrezId
	a.Transcript = trans.Transcript
}

func (a *GeneAnnoItem) SetRegion(region string) {
	a.Region = region
}

func (a *GeneAnnoItem) SetEvent(event string) {
	a.Event = event
}

func (a *GeneAnnoItem) SetNAChange(change string) {
	a.NAChange = change
}

func (a *GeneAnnoItem) SetAAChange(change string) {
	a.AAChange = change
}

func (a *GeneAnnoItem) AnnoInGene(snp variant.Snv, trans transcript.Transcript) {
	a.GeneSymbol = trans.Gene
	a.GeneEntrezId = trans.EntrezId
	a.Transcript = trans.Transcript
	regionIndex, exonLen := trans.Regions.FindOne(snp.Start, trans.Strand)
	region := trans.Regions[regionIndex]
	if region.Type == "intron" {
		a.AnnoInIntron(snp, trans, regionIndex, exonLen)
	} else if strings.HasPrefix(region.Type, "UTR") {
		a.AnnoInUTR(trans, regionIndex)
	} else {
		a.AnnoInCDS(snp, trans, regionIndex, exonLen)
	}
}
