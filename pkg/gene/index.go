package gene

import (
	"bufio"
	"open-anno/pkg"
	"os"
	"sort"
)

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

func (this TransIndexes) Len() int           { return len(this) }
func (this TransIndexes) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this TransIndexes) Less(i, j int) bool { return this[i].Start < this[j].Start }

func NewTransIndexes(step int) TransIndexes {
	indexes := make(TransIndexes, 0)
	for chrom, length := range GENOME {
		for start := 1; start <= length; start += step {
			end := start + step - 1
			if end > length {
				end = length
			}
			index := TransIndex{Chrom: chrom, Start: start, End: end}
			indexes = append(indexes, index)
		}
	}
	sort.Sort(indexes)
	return indexes
}

func ReadTransIndexDB(dbFile string) (TransIndexes, error) {
	indexes := make(TransIndexes, 0)
	fi, err := os.Open(dbFile)
	if err != nil {
		return indexes, err
	}
	defer fi.Close()
	scanner := bufio.NewScanner(fi)
	for scanner.Scan() {
		index := pkg.FromJSON[TransIndex](scanner.Text())
		indexes = append(indexes, index)
	}
	return indexes, err
}
