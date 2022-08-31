package io

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/scheme"
	"strconv"
	"strings"
)

type DBVarScanner struct {
	Scanner[scheme.DBVar]
	Header string
}

func NewDBVarScanner(reader Reader) DBVarScanner {
	scanner := NewScanner[scheme.DBVar](reader)
	scanner.Scan()
	text := scanner.Text()
	header := strings.TrimLeft(text, "#")
	return DBVarScanner{Scanner: scanner, Header: header}
}

func (this DBVarScanner) Row() (scheme.DBVar, error) {
	text := this.Text()
	fields := strings.Split(text, "\t")
	dbvar := scheme.DBVar{
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
	Scanner[scheme.DBVarIdx]
	BinSize int
}

func (this DBVarIdxScanner) Row() (scheme.DBVarIdx, error) {
	field := strings.Split(this.Text(), "\t")
	bin := fmt.Sprintf("%s\t%s", field[0], field[1])
	start, err := strconv.ParseInt(field[2], 10, 0)
	if err != nil {
		return scheme.DBVarIdx{}, err
	}
	end, err := strconv.ParseInt(field[3], 10, 0)
	if err != nil {
		return scheme.DBVarIdx{}, err
	}
	return scheme.DBVarIdx{Bin: bin, Start: start, End: end}, err
}

func NewDBVarIdxScanner(reader Reader) DBVarIdxScanner {
	scanner := NewScanner[scheme.DBVarIdx](reader)
	scanner.Scan()
	binSize, err := strconv.Atoi(strings.Split(scanner.Text(), "\t")[1])
	if err != nil {
		log.Fatal(err)
	}
	return DBVarIdxScanner{Scanner: scanner, BinSize: binSize}
}
