package io

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/schema"
	"strconv"
	"strings"
)

type DBVarScanner struct {
	Scanner[schema.DBVar]
	Header string
}

func NewDBVarScanner(reader Reader) DBVarScanner {
	scanner := NewScanner[schema.DBVar](reader)
	scanner.Scan()
	text := scanner.Text()
	header := strings.TrimLeft(text, "#")
	return DBVarScanner{Scanner: scanner, Header: header}
}

func (this DBVarScanner) Row() (schema.DBVar, error) {
	text := this.Text()
	fields := strings.Split(text, "\t")
	dbvar := schema.DBVar{
		Chrom: pkg.FormatChrom(fields[0]),
		Ref:   fields[3],
		Alt:   fields[4],
		Text:  text,
	}
	var err error
	dbvar.Start, err = strconv.Atoi(fields[1])
	if err != nil {
		return dbvar, err
	}
	dbvar.End, err = strconv.Atoi(fields[2])
	if err != nil {
		return dbvar, err
	}
	return dbvar, nil
}

type DBVarIdxScanner struct {
	Scanner[schema.DBVarIdx]
	BinSize int
}

func (this DBVarIdxScanner) Row() (schema.DBVarIdx, error) {
	field := strings.Split(this.Text(), "\t")
	bin := fmt.Sprintf("%s\t%s", field[0], field[1])
	start, err := strconv.ParseInt(field[2], 10, 0)
	if err != nil {
		return schema.DBVarIdx{}, err
	}
	end, err := strconv.ParseInt(field[3], 10, 0)
	if err != nil {
		return schema.DBVarIdx{}, err
	}
	return schema.DBVarIdx{Bin: bin, Start: start, End: end}, err
}

func NewDBVarIdxScanner(reader Reader) DBVarIdxScanner {
	scanner := NewScanner[schema.DBVarIdx](reader)
	scanner.Scan()
	binSize, err := strconv.Atoi(strings.Split(scanner.Text(), "\t")[1])
	if err != nil {
		log.Fatal(err)
	}
	return DBVarIdxScanner{Scanner: scanner, BinSize: binSize}
}
