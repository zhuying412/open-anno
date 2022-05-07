package io

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"open-anno/pkg"
	"os"
	"strconv"
	"strings"
)

type BED struct {
	Chrom string `json:"chrom"`
	Start int    `json:"start"`
	End   int    `json:"end"`
	Name  string `json:"name"`
}

type BEDs []BED

func (this BEDs) Len() int      { return len(this) }
func (this BEDs) Swap(i, j int) { this[i], this[j] = this[j], this[i] }
func (this BEDs) Less(i, j int) bool {
	if this[i].Start == this[j].Start {
		return this[i].End < this[j].End
	}
	return this[i].Start < this[j].Start
}

type BEDScanner struct {
	scanner *bufio.Scanner
}

func NewBEDScanner(reader io.Reader) BEDScanner {
	scanner := bufio.NewScanner(reader)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)
	return BEDScanner{scanner: scanner}
}

func (this *BEDScanner) Scan() bool {
	return this.scanner.Scan()
}

func (this BEDScanner) Text() string {
	return this.scanner.Text()
}

func (this BEDScanner) Row() (BED, error) {
	fields := strings.Split(this.Text(), "\t")
	if len(fields) > 4 {
		log.Printf("Error: the column is %d, %d-%d will be ignored\n", len(fields), 4, len(fields))
	}
	bed := BED{
		Chrom: pkg.FormatChrom(fields[0]),
		Name:  fields[3],
	}
	var err error
	bed.Start, err = strconv.Atoi(fields[1])
	if err != nil {
		return bed, err
	}
	bed.End, err = strconv.Atoi(fields[2])
	if err != nil {
		return bed, err
	}
	return bed, nil
}

func ReadBEDs(infile string) (BEDs, error) {
	var beds BEDs
	reader, err := os.Open(infile)
	if err != nil {
		return beds, err
	}
	defer reader.Close()
	scanner := NewBEDScanner(reader)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return beds, err
		}
		beds = append(beds, row)
	}
	return beds, err
}

func WriteBEDs(outfile string, beds ...BEDs) error {
	writer, err := os.Create(outfile)
	if err != nil {
		return err
	}
	for _, rows := range beds {
		for _, row := range rows {
			fmt.Fprintf(writer, "%s\t%d\t%d\t%s\n", row.Chrom, row.Start, row.End, row.Name)
		}
	}
	return nil
}
