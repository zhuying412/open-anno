package io

import (
	"bufio"
	"io"
	"os"
	"strings"
)

type DBVarScanner struct {
	VarScanner
	Header string
}

func NewDBVarScanner(reader io.Reader) DBVarScanner {
	scanner := bufio.NewScanner(reader)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)
	scanner.Scan()
	header := strings.TrimLeft(scanner.Text(), "#")
	return DBVarScanner{VarScanner: VarScanner{scanner: scanner}, Header: header}
}

type DBBEDScanner struct {
	BEDScanner
	Header string
}

func NewDBBEDScanner(reader io.Reader) DBBEDScanner {
	scanner := bufio.NewScanner(reader)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)
	scanner.Scan()
	header := strings.TrimLeft(scanner.Text(), "#")
	return DBBEDScanner{BEDScanner: BEDScanner{scanner: scanner}, Header: header}
}

func ReadDBBEDs(infile string) (BEDs, string, error) {
	var beds BEDs
	reader, err := os.Open(infile)
	if err != nil {
		return beds, "", err
	}
	defer reader.Close()
	scanner := NewDBBEDScanner(reader)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return beds, scanner.Header, err
		}
		beds = append(beds, row)
	}
	return beds, scanner.Header, err
}
