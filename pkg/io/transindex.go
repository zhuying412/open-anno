package io

import (
	"fmt"
	"open-anno/pkg/scheme"
	"open-anno/pkg/seq"
	"sort"
	"strings"
)

func CreateTransIndexes(transcripts scheme.Transcripts, step int, outfile string) error {
	transIndexes := make(scheme.TransIndexes, 0)
	for chrom, length := range seq.GENOME {
		for start := 1; start <= length; start += step {
			end := start + step - 1
			if end > length {
				end = length
			}
			transIndex := scheme.TransIndex{Chrom: chrom, Start: start, End: end}
			transIndexes = append(transIndexes, transIndex)
		}
	}
	sort.Sort(transIndexes)
	writer, err := NewIoWriter(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	for _, transIndex := range transIndexes {
		transIndex.SetTranscripts(transcripts)
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

func ReadTransIndexs(dbFile string) (scheme.TransIndexes, error) {
	indexes := make(scheme.TransIndexes, 0)
	reader, err := NewIoReader(dbFile)
	if err != nil {
		return indexes, err
	}
	defer reader.Close()
	scanner := NewBEDScanner(reader)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return indexes, err
		}
		indexes = append(indexes, scheme.TransIndex{
			Chrom:       row.Chrom,
			Start:       row.Start,
			End:         row.End,
			Transcripts: strings.Split(row.Name, ","),
		})
	}
	sort.Sort(indexes)
	return indexes, err
}