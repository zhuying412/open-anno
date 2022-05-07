package io

import (
	"bytes"
	"fmt"
	"log"
	"os"
	"strings"
)

type DBAnno struct {
	Length int
	Header string
	Data   map[string]string
}

func (this DBAnno) Text(id string) string {
	if text, ok := this.Data[id]; ok {
		return text
	}
	var buffer bytes.Buffer
	for i := 0; i < this.Length; i++ {
		if i == 0 {
			buffer.WriteString(".")
		} else {
			buffer.WriteString("\t.")
		}
	}
	return buffer.String()
}

type DBAnnos []DBAnno

func (this DBAnnos) Header() string {
	var buffer bytes.Buffer
	for i, dbAnno := range this {
		if i == 0 {
			buffer.WriteString(dbAnno.Header)
		} else {
			buffer.WriteString("\t" + dbAnno.Header)
		}
	}
	return buffer.String()
}

func (this DBAnnos) Text(id string) string {
	var buffer bytes.Buffer
	for i, dbAnno := range this {
		if i == 0 {
			buffer.WriteString(dbAnno.Text(id))
		} else {
			buffer.WriteString("\t" + dbAnno.Text(id))
		}
	}
	return buffer.String()
}

func readOtherbased(infile string) (DBAnno, error) {
	var dbAnno DBAnno
	reader, err := os.Open(infile)
	if err != nil {
		return dbAnno, err
	}
	defer reader.Close()
	scanner := NewDBVarScanner(reader)
	headers := strings.Split(scanner.Header, "\t")[5:]
	dbAnno = DBAnno{
		Header: strings.Join(headers, "\t"),
		Length: len(headers),
		Data:   map[string]string{},
	}
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return dbAnno, err
		}
		dbAnno.Data[row.ID()] = strings.Join(row.Otherinfo, "\t")
	}
	return dbAnno, err
}

func readAVinput(infile string) (DBAnno, error) {
	var dbAnno DBAnno
	reader, err := os.Open(infile)
	if err != nil {
		return dbAnno, err
	}
	defer reader.Close()
	scanner := NewVarScanner(reader)
	flag := false
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return dbAnno, err
		}
		if !flag {
			headers := make([]string, len(row.Otherinfo))
			for i := 0; i < len(row.Otherinfo); i++ {
				headers[i] = fmt.Sprintf("Otherinfo%d", i+1)
			}
			dbAnno.Header = strings.Join(headers, "\t")
			dbAnno.Length = len(row.Otherinfo)
			dbAnno.Data = map[string]string{}
			flag = true
		}
		dbAnno.Data[row.ID()] = strings.Join(row.Otherinfo, "\t")
	}
	return dbAnno, err
}

func MergeAnno(outfile string, avinput string, genebased string, otherbaseds ...string) error {
	dbAnnos := make(DBAnnos, 0)
	for _, otherbased := range otherbaseds {
		dbAnno, err := readOtherbased(otherbased)
		if err != nil {
			return err
		}
		dbAnnos = append(dbAnnos, dbAnno)
	}
	dbAnno, err := readAVinput(avinput)
	if err != nil {
		return err
	}
	dbAnnos = append(dbAnnos, dbAnno)
	reader, err := os.Open(genebased)
	if err != nil {
		log.Fatal(err)
	}
	defer reader.Close()
	writer, err := os.Create(outfile)
	if err != nil {
		log.Fatal(err)
	}
	defer writer.Close()
	scanner := NewDBVarScanner(reader)
	fmt.Fprintf(writer, "%s\t%s\n", scanner.Header, dbAnnos.Header())
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return err
		}
		fmt.Fprintf(writer, "%s\t%s\n", scanner.Text(), dbAnnos.Text(row.ID()))
	}
	return nil
}
