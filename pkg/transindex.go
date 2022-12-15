package pkg

import (
	"fmt"
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

func CreateTransIndexes(gpes GenePreds, genome map[string]int, step int, idxFile string) error {
	transIndexes := make(TransIndexes, 0)
	for chrom, length := range genome {
		for start := 1; start <= length; start += step {
			end := start + step - 1
			if end > length {
				end = length
			}
			transIndex := TransIndex{Chrom: chrom, Start: start, End: end}
			transIndexes = append(transIndexes, transIndex)
		}
	}
	sort.Sort(transIndexes)
	writer, err := NewIOWriter(idxFile)
	if err != nil {
		return err
	}
	defer writer.Close()
	for _, transIndex := range transIndexes {
		for _, gpe := range gpes {
			if transIndex.Chrom == gpe.Chrom && transIndex.Start <= gpe.TxEnd && transIndex.End >= gpe.TxStart {
				transIndex.Transcripts = append(transIndex.Transcripts, gpe.PK())
			}
		}
		if len(transIndex.Transcripts) > 0 {
			fmt.Fprintf(writer, "%s\t%d\t%d\t%s\n",
				transIndex.Chrom,
				transIndex.Start,
				transIndex.End,
				strings.Join(transIndex.Transcripts, ","),
			)
		}
	}
	return err
}

func ReadTransIndexes(dbFile string) (TransIndexes, error) {
	indexes := make(TransIndexes, 0)
	reader, err := NewIOReader(dbFile)
	if err != nil {
		return indexes, err
	}
	defer reader.Close()
	scanner := NewIOScanner(reader)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		index := TransIndex{
			Chrom:       row[0],
			Transcripts: strings.Split(row[3], ","),
		}
		index.Start, err = strconv.Atoi(row[1])
		if err != nil {
			return indexes, err
		}
		index.End, err = strconv.Atoi(row[2])
		if err != nil {
			return indexes, err
		}
		indexes = append(indexes, index)
	}
	sort.Sort(indexes)
	return indexes, err
}
