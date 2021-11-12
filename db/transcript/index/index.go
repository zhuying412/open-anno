package index

type TranscriptIndex struct {
	Chrom       string   `json:"chrom"`
	Start       int      `json:"start"`
	End         int      `json:"end"`
	Transcripts []string `json:"transcripts"`
}

func (t *TranscriptIndex) SetTranscript(transcripts Transcripts) {
	for _, trans := range transcripts {
		if t.Chrom == trans.Chrom && t.Start <= trans.End && t.End >= trans.Start {
			t.Transcripts = append(t.Transcripts, trans.SN())
		}
	}
}

type TranscriptIndexes []TranscriptIndex

func (t TranscriptIndexes) Len() int {
	return len(t)
}

func (t TranscriptIndexes) Less(i, j int) bool {
	return t[i].Start < t[j].Start || t[i].End < t[j].End
}

func (t TranscriptIndexes) Swap(i, j int) {
	t[i], t[j] = t[j], t[i]
}

func (t TranscriptIndexes) FilterByChrom(chrom string) TranscriptIndexes {
	transIndexes := make(TranscriptIndexes, 0)
	for _, transIndex := range t {
		if transIndex.Chrom == chrom {
			transIndexes = append(transIndexes, transIndex)
		}
	}
	return transIndexes
}
