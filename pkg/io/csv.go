package io

import (
	"bufio"
	"io"
	"strings"
)

type CSVScanner struct {
	scanner    *bufio.Scanner
	FieldNames []string
}

func NewCSVScanner(reader io.Reader) CSVScanner {
	scanner := bufio.NewScanner(reader)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)
	scanner.Scan()
	return CSVScanner{scanner: scanner, FieldNames: strings.Split(scanner.Text(), "\t")}
}

func (this *CSVScanner) Scan() bool {
	return this.scanner.Scan()
}

func (this CSVScanner) Text() string {
	return this.scanner.Text()
}

func (this CSVScanner) Row() map[string]string {
	row := make(map[string]string)
	fields := strings.Split(this.Text(), "\t")
	for i, header := range this.FieldNames {
		row[header] = fields[i]
	}
	return row
}
