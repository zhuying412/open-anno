package scheme

import (
	"bytes"
	"fmt"
	"open-anno/pkg/seq"

	"github.com/brentp/faidx"
)

// Transcript Refgene
type Transcript struct {
	Name       string  `json:"name"`
	Chrom      string  `json:"chrom"`
	Strand     string  `json:"strand"`
	TxStart    int     `json:"txStart"`
	TxEnd      int     `json:"txEnd"`
	CdsStart   int     `json:"cdsStart"`
	CdsEnd     int     `json:"cdsEnd"`
	ExonCount  int     `json:"exonCount"`
	ExonStarts []int   `json:"exonStarts"`
	ExonEnds   []int   `json:"exonEnds"`
	Gene       string  `json:"gene"`
	GeneID     string  `json:"gene_id"`
	CdsStat    string  `json:"cdsStartStat"`
	Regions    Regions `json:"regions"`
}

func (this Transcript) CdsCount() int {
	var count int
	for _, region := range this.Regions {
		if region.Type == RType_CDS {
			count++
		}
	}
	return count
}

func (this Transcript) CLen() int {
	var length int
	for _, region := range this.Regions {
		if region.Type == RType_CDS {
			length += region.End - region.Start + 1
		}
	}
	return length
}

func (this Transcript) CDNA() string {
	var cdna bytes.Buffer
	for _, region := range this.Regions {
		if region.Type == RType_CDS {
			cdna.WriteString(region.Sequence)
		}
	}
	return cdna.String()
}

func (this Transcript) DNA() string {
	var cdna bytes.Buffer
	for _, region := range this.Regions {
		cdna.WriteString(region.Sequence)
	}
	return cdna.String()
}

func (this Transcript) IsUnk() bool {
	return this.CdsEnd-this.CdsStart+1 == 0
}

func (this *Transcript) SetRegions(mrna *faidx.Faidx, geneInfoMap GeneInfoMap, seqRequired bool) error {
	regions, err := NewRegions(*this)
	if err != nil {
		return err
	}
	for i, region := range regions {
		mrnaName := fmt.Sprintf("%s:%s:%s", this.Chrom, this.Gene, this.Name)
		regions[i].Sequence, err = seq.Fetch(mrna, mrnaName, region.Start-this.TxStart, region.End-this.TxStart+1)
		if err != nil && seqRequired {
			return err
		}
	}
	this.Regions = regions
	if geneInfo, ok := geneInfoMap[this.Chrom][this.Gene]; ok {
		this.GeneID = geneInfo.EntrezId
	}
	return nil
}

// Transcripts
type Transcripts map[string]Transcript

func (this Transcripts) FilterChrom(chrom string, geneInfoMap GeneInfoMap, mrna *faidx.Faidx, seqRequired bool) (Transcripts, error) {
	transcripts := make(Transcripts)
	for sn, trans := range this {
		if trans.Chrom == chrom {
			err := trans.SetRegions(mrna, geneInfoMap, seqRequired)
			if err != nil {
				return transcripts, err
			}
			transcripts[sn] = trans
		}
	}
	return transcripts, nil
}
