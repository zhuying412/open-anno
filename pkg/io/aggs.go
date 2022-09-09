package io

import (
	"fmt"
	"open-anno/pkg/schema"
	"os"
	"sort"
	"strings"
)

type SnvTransAnno struct {
	Snv    string
	Gene   string
	GeneID string
	Region string
	Event  string
	Detail string
}

func (this SnvTransAnno) ID() string {
	return fmt.Sprintf("%s\t%s", this.Snv, this.Gene)
}

type SnvGBAnnoScanner struct {
	GenericsTSVScanner[SnvTransAnno]
	DBname string
}

func NewSnvGBAnnoScanner(reader Reader, dbname string) SnvGBAnnoScanner {
	scanner := NewGenericsTSVScanner[SnvTransAnno](reader)
	return SnvGBAnnoScanner{GenericsTSVScanner: scanner, DBname: dbname}
}

func (this SnvGBAnnoScanner) Row() (SnvTransAnno, error) {
	row, err := this.StringMapRow()
	if err != nil {
		return SnvTransAnno{}, err
	}
	return SnvTransAnno{
		Snv:    fmt.Sprintf("%s\t%s\t%s\t%s\t%s", row["Chr"], row["Start"], row["End"], row["Ref"], row["Alt"]),
		Gene:   row[this.DBname+".Gene"],
		GeneID: row[this.DBname+".GeneID"],
		Event:  row[this.DBname+".Event"],
		Region: row[this.DBname+".Region"],
		Detail: row[this.DBname+".Detail"],
	}, nil
}

func AggsSnvGBAnno(annoGBOutput, dbname, outfile string) error {
	reader, err := NewIoReader(annoGBOutput)
	if err != nil {
		return err
	}
	defer reader.Close()
	scanner := NewSnvGBAnnoScanner(reader, dbname)
	snvGeneAnnos := make(map[string]schema.SnvGeneAnno)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return err
		}
		snvGeneAnno, ok := snvGeneAnnos[row.ID()]
		if !ok {
			snvGeneAnno = schema.SnvGeneAnno{
				Snv:    row.Snv,
				Gene:   row.Gene,
				GeneID: row.GeneID,
			}
		}
		snvGeneAnno.AddInfo(row.Event, row.Region, row.Detail)
		snvGeneAnnos[row.ID()] = snvGeneAnno
	}
	snvGeneAnnoTexts := make([]string, len(snvGeneAnnos))
	i := 0
	for _, snvGeneAnno := range snvGeneAnnos {
		snvGeneAnnoTexts[i] = snvGeneAnno.Text()
		i++
	}
	sort.Strings(snvGeneAnnoTexts)
	writer, err := os.Create(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprintf(writer, "%s\n", strings.Join(scanner.FieldNames, "\t"))
	for _, text := range snvGeneAnnoTexts {
		fmt.Fprintf(writer, "%s\n", text)
	}
	return nil
}
