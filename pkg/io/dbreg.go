package io

import (
	"open-anno/pkg"
	"open-anno/pkg/scheme"
	"strconv"
	"strings"
)

type DBRegScanner struct {
	Scanner[scheme.DBReg]
	FieldName string
}

func NewDBRegScanner(reader Reader) DBRegScanner {
	scanner := NewScanner[scheme.DBReg](reader)
	scanner.Scan()
	return DBRegScanner{Scanner: scanner, FieldName: strings.Split(scanner.Text(), "\t")[3]}
}

func (this DBRegScanner) Row() (scheme.DBReg, error) {
	fields := strings.Split(this.Text(), "\t")
	dbreg := scheme.DBReg{
		Chrom: pkg.FormatChrom(fields[0]),
		Info:  fields[3],
	}
	var err error
	dbreg.Start, err = strconv.Atoi(fields[1])
	if err != nil {
		return dbreg, err
	}
	dbreg.End, err = strconv.Atoi(fields[2])
	if err != nil {
		return dbreg, err
	}
	return dbreg, nil
}
