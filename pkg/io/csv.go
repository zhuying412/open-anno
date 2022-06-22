package io

import (
	"errors"
	"fmt"
	"io"
	"strings"
)

type CSVScanner struct {
	Scanner[map[string]string]
	FieldNames []string
}

func NewCSVScanner(reader io.Reader) CSVScanner {
	scanner := NewScanner[map[string]string](reader)
	scanner.Scan()
	return CSVScanner{Scanner: scanner, FieldNames: strings.Split(scanner.Text(), "\t")}
}

func (this CSVScanner) Row() (map[string]string, error) {
	row := make(map[string]string)
	fields := strings.Split(this.Text(), "\t")
	if len(fields) == len(this.FieldNames) {
		return row, errors.New(fmt.Sprintf("the column is not equal the FieldNames: %s", this.Text()))
	}
	for i, header := range this.FieldNames {
		row[header] = fields[i]
	}
	return row, nil
}
