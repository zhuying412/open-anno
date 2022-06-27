package io

import (
	"fmt"
	"io"
	"log"
	"strconv"
	"strings"
)

type DBVarIdx struct {
	Bin   string `json:"bin"`
	Start int64  `json:"start"`
	End   int64  `json:"end"`
}
type DBVarIdxScanner struct {
	Scanner[DBVarIdx]
	BinSize int
}

func (this DBVarIdxScanner) Row() (DBVarIdx, error) {
	field := strings.Split(this.Text(), "\t")
	bin := fmt.Sprintf("%s:%s", field[0], field[1])
	start, err := strconv.ParseInt(field[2], 10, 0)
	if err != nil {
		return DBVarIdx{}, err
	}
	end, err := strconv.ParseInt(field[3], 10, 0)
	if err != nil {
		return DBVarIdx{}, err
	}
	return DBVarIdx{Bin: bin, Start: start, End: end}, err
}

func NewDBVarIdxScanner(reader io.Reader) DBVarIdxScanner {
	scanner := NewScanner[DBVarIdx](reader)
	scanner.Scan()
	binSize, err := strconv.Atoi(strings.Split(scanner.Text(), "\t")[1])
	if err != nil {
		log.Fatal(err)
	}
	return DBVarIdxScanner{Scanner: scanner, BinSize: binSize}
}

type DBVarScanner struct {
	VarScanner
	Header string
}

func NewDBVarScanner(reader io.Reader) DBVarScanner {
	scanner := NewVarScanner(reader)
	scanner.Scan()
	text := scanner.Text()
	header := strings.TrimLeft(text, "#")
	return DBVarScanner{VarScanner: scanner, Header: header}
}

type DBRegScanner struct {
	BEDScanner
	Header string
}

func NewDBRegScanner(reader io.Reader) DBRegScanner {
	scanner := NewBEDScanner(reader)
	scanner.Scan()
	text := scanner.Text()
	header := strings.TrimLeft(text, "#")
	return DBRegScanner{BEDScanner: scanner, Header: header}
}
