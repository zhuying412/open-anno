package gene

import (
	"bufio"
	"os"
	"sort"
	"strconv"
	"strings"
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

func ReadTransIndexs(dbFile string) (TransIndexes, error) {
	indexes := make(TransIndexes, 0)
	fi, err := os.Open(dbFile)
	if err != nil {
		return indexes, err
	}
	defer fi.Close()
	scanner := bufio.NewScanner(fi)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return indexes, err
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return indexes, err
		}
		indexes = append(indexes, TransIndex{
			Chrom:       fields[0],
			Start:       start,
			End:         end,
			Transcripts: strings.Split(fields[3], ","),
		})
	}
	sort.Sort(indexes)
	return indexes, err
}
