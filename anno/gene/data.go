package gene

import (
	"log"
	"open-anno/pkg/io"
	"open-anno/pkg/scheme"

	"github.com/brentp/faidx"
)

type GeneData struct {
	Transcripts  scheme.Transcripts                    `json:"transcirpts"`
	TransIndexes scheme.TransIndexes                   `json:"transindexes"`
	GeneInfo     map[string]map[string]scheme.GeneInfo `json:"symboltoid"`
	MrnaFaidx    *faidx.Faidx                          `json:"mrna"`
}

func (this GeneData) FilterTranscripts(chrom string, seqRequired bool) (scheme.Transcripts, error) {
	return this.Transcripts.FilterChrom(chrom, this.GeneInfo, this.MrnaFaidx, seqRequired)
}

func (this GeneData) FilterTransIndexes(chrom string) scheme.TransIndexes {
	return this.TransIndexes.FilterChrom(chrom)
}

func NewGeneData(refgeneFile, indexFile, geneInfoFile, mrnaFile string) (GeneData, error) {
	var geneData GeneData
	var err error
	log.Printf("Read GenePred: %s ...", refgeneFile)
	geneData.Transcripts, err = io.ReadGenePred(refgeneFile)
	if err != nil {
		return geneData, err
	}
	log.Printf("Read GeneInfo from %s", geneInfoFile)
	geneData.GeneInfo, err = io.ReadGeneInfo(geneInfoFile)
	if err != nil {
		return geneData, err
	}
	log.Printf("Read GenePred Index: %s ...", indexFile)
	geneData.TransIndexes, err = io.ReadTransIndexs(indexFile)
	if err != nil {
		return geneData, err
	}
	if mrnaFile != "" {
		log.Printf("Read mRNA: %s ...", mrnaFile)
		geneData.MrnaFaidx, err = faidx.New(mrnaFile)
		if err != nil {
			return geneData, err
		}
	}
	return geneData, err
}
