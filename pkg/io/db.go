package io

import (
	"io"
	"strings"
)

type DBVarScanner struct {
	VarScanner
	Header string
}

func NewDBVarScanner(reader io.Reader) DBVarScanner {
	scanner := NewVarScanner(reader)
	scanner.Scan()
	header := strings.TrimLeft(scanner.Text(), "#")
	return DBVarScanner{VarScanner: scanner, Header: header}
}

type DBRegScanner struct {
	BEDScanner
	Header string
}

func NewDBRegScanner(reader io.Reader) DBRegScanner {
	scanner := NewBEDScanner(reader)
	scanner.Scan()
	header := strings.TrimLeft(scanner.Text(), "#")
	return DBRegScanner{BEDScanner: scanner, Header: header}
}

func ReadDBRegs(infile string) (BEDs, string, error) {
	var beds BEDs
	reader, err := NewIoReader(infile)
	if err != nil {
		return beds, "", err
	}
	defer reader.Close()
	var scanner DBRegScanner
	scanner = NewDBRegScanner(reader)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return beds, scanner.Header, err
		}
		beds = append(beds, row)
	}
	return beds, scanner.Header, err
}
