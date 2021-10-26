package db

import (
	"grandanno/seq"
)

func NewMrnaFastaMap(refgenes Refgenes, chromSeq seq.Fasta) map[string]seq.Fasta {
	mrna := make(map[string]seq.Fasta)
	for _, refgene := range refgenes {
		if sequence, ok := chromSeq[refgene.Chrom]; ok {
			if _, ok = mrna[refgene.Chrom]; !ok {
				mrna[refgene.Chrom] = make(seq.Fasta)
			}
			mrna[refgene.Chrom][refgene.SN()] = sequence.SubSeq(refgene.Start-1, refgene.End-refgene.Start+1)
		}
	}
	return mrna
}
