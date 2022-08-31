package scheme

type TransIndex struct {
	Chrom       string   `json:"chrom"`
	Start       int      `json:"start"`
	End         int      `json:"end"`
	Transcripts []string `json:"transcript"`
}

func (this *TransIndex) SetTranscripts(transcripts Transcripts) {
	for _, trans := range transcripts {
		if this.Chrom == trans.Chrom && this.Start <= trans.TxEnd && this.End >= trans.TxStart {
			this.Transcripts = append(this.Transcripts, trans.Name)
		}
	}
}

type TransIndexes []TransIndex

func (this TransIndexes) Len() int      { return len(this) }
func (this TransIndexes) Swap(i, j int) { this[i], this[j] = this[j], this[i] }
func (this TransIndexes) Less(i, j int) bool {
	if this[i].Chrom != this[j].Chrom {
		return this[i].Chrom < this[j].Chrom
	}
	return this[i].Start < this[j].Start
}

func (this TransIndexes) FilterChrom(chrom string) TransIndexes {
	indexes := make(TransIndexes, 0)
	for _, index := range this {
		if index.Chrom == chrom {
			indexes = append(indexes, index)
		}
	}
	return indexes
}
