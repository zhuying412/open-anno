package gene

import (
	"log"
	"open-anno/pkg/io"
	"open-anno/pkg/io/refgene"

	"github.com/brentp/faidx"
)

type GeneData struct {
	Transcripts  refgene.Transcripts  `json:"transcirpts"`
	TransIndexes refgene.TransIndexes `json:"transindexes"`
	SymbolToId   map[string]string    `json:"symboltoid"`
	MrnaFaidx    *faidx.Faidx         `json:"mrna"`
}

func (this GeneData) FilterTranscripts(chrom string, seqRequired bool) (refgene.Transcripts, error) {
	return this.Transcripts.FilterChrom(chrom, this.SymbolToId, this.MrnaFaidx, seqRequired)
}

func (this GeneData) FilterTransIndexes(chrom string) refgene.TransIndexes {
	return this.TransIndexes.FilterChrom(chrom)
}

func NewGeneData(refgeneFile, indexFile, ncbiGeneInfoFile, gene2RefseqFile, mrnaFile string) (GeneData, error) {
	var geneData GeneData
	var err error
	log.Printf("Read GenePred: %s ...", refgeneFile)
	geneData.Transcripts, err = refgene.ReadRefgene(refgeneFile)
	if err != nil {
		return geneData, err
	}
	log.Printf("New Gene Symbol to EntrezID from %s, %s", ncbiGeneInfoFile, gene2RefseqFile)
	geneData.SymbolToId, err = io.NewGeneSymbolToId(gene2RefseqFile, ncbiGeneInfoFile, refgeneFile)
	if err != nil {
		return geneData, err
	}
	log.Printf("Read GenePred Index: %s ...", indexFile)
	geneData.TransIndexes, err = refgene.ReadTransIndexs(indexFile)
	if err != nil {
		return geneData, err
	}
	if mrnaFile == "" {
		log.Printf("Read mRNA: %s ...", mrnaFile)
		geneData.MrnaFaidx, err = faidx.New(mrnaFile)
		if err != nil {
			return geneData, err
		}
	}
	return geneData, err
}
