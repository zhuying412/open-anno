package io

import (
	"errors"
	"fmt"
	"log"
	"strings"
)

type GenericsTSVScanner[T any] struct {
	Scanner[T]
	FieldNames []string
}

func NewGenericsTSVScanner[T any](reader Reader) GenericsTSVScanner[T] {
	scanner := NewScanner[T](reader)
	scanner.Scan()
	return GenericsTSVScanner[T]{Scanner: scanner, FieldNames: strings.Split(scanner.Text(), "\t")}
}

func (this GenericsTSVScanner[T]) StringMapRow() (map[string]string, error) {
	row := make(map[string]string)
	fields := strings.Split(this.Text(), "\t")
	if len(fields) != len(this.FieldNames) {
		return row, errors.New(fmt.Sprintf("the column is not equal the FieldNames: %s", this.Text()))
	}
	for i, header := range this.FieldNames {
		row[header] = fields[i]
	}
	return row, nil
}

func (this GenericsTSVScanner[T]) Row() (T, error) {
	log.Fatal("Not implemented")
	return *new(T), nil
}

type CSVScanner struct {
	GenericsTSVScanner[map[string]string]
}

func NewCSVScanner(reader Reader) CSVScanner {
	scanner := NewGenericsTSVScanner[map[string]string](reader)
	return CSVScanner{GenericsTSVScanner: scanner}
}

func (this CSVScanner) Row() (map[string]string, error) {
	return this.StringMapRow()
}
