package index

import (
	"fmt"
)

type Transcript struct {
	ID         string `json:"id"`
	Chrom      string `json:"chrom"`
	Start      int    `json:"start"`
	End        int    `json:"end"`
	UpStream   int    `json:"up_stream"`
	DownStream int    `json:"down_stream"`
}

func (t Transcript) SN() string {
	return fmt.Sprintf("%s|%s:%d:%d", t.ID, t.Chrom, t.Start, t.End)
}

type Transcripts []Transcript

func (Transcripts Transcripts) Len() int {
	return len(Transcripts)
}

func (t Transcripts) Less(i, j int) bool {
	return t[i].Start < t[j].Start || t[i].End < t[j].End
}

func (t Transcripts) Swap(i, j int) {
	t[i], t[j] = t[j], t[i]
}

func (t Transcripts) FilterByChrom(chrom string) Transcripts {
	trans := make(Transcripts, 0)
	for _, transIndex := range t {
		if transIndex.Chrom == chrom {
			trans = append(trans, transIndex)
		}
	}
	return trans
}
