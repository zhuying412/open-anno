package snv

import (
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"os"
	"sort"
	"strings"
)

type GeneAnno struct {
	Snv     string
	Gene    string
	GeneID  string
	Regions []string
	Events  []string
	Details []string
}

func (this *GeneAnno) AddInfo(event, region, detail string) {
	if event != "." && pkg.FindArr(this.Events, event) < 0 {
		this.Events = append(this.Events, event)
	}
	if region != "." && pkg.FindArr(this.Regions, region) < 0 {
		this.Regions = append(this.Regions, region)
	}
	if detail != "." && pkg.FindArr(this.Details, detail) < 0 {
		this.Details = append(this.Details, detail)
	}
}

func (this GeneAnno) Region() string {
	var regions1, regions2, regions3, regions4 []string
	for _, region := range this.Regions {
		switch region {
		case "exonic", "splicing", "exonic_splicing", "transcript":
			regions1 = append(regions1, region)
		case "ncRNA":
			regions2 = append(regions2, region)
		case "UTR3", "UTR5":
			regions3 = append(regions3, region)
		case "intronic":
			regions4 = append(regions4, region)
		}
	}
	if len(regions1) > 0 {
		return strings.Join(regions1, ",")
	}
	if len(regions2) > 0 {
		return strings.Join(regions2, ",")
	}
	if len(regions3) > 0 {
		return strings.Join(regions3, ",")
	}
	if len(regions4) > 0 {
		return strings.Join(regions4, ",")
	}
	return "."
}

func (this GeneAnno) Event() string {
	if len(this.Events) > 0 {
		return strings.Join(this.Events, ",")
	}
	return "."
}

func (this GeneAnno) Detail() string {
	if len(this.Details) > 0 {
		return strings.Join(this.Details, ",")
	}
	return "."
}

func (this GeneAnno) Text() string {
	return fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s", this.Snv, this.Gene, this.GeneID, this.Event(), this.Region(), this.Detail())
}

type TransAnnoRow struct {
	Snv    string
	Gene   string
	GeneID string
	Region string
	Event  string
	Detail string
}

func (this TransAnnoRow) ID() string {
	return fmt.Sprintf("%s\t%s", this.Snv, this.Gene)
}

type TransAnnoScanner struct {
	io.GenericsTSVScanner[TransAnno]
	DBname string
}

func NewTransAnnoScanner(reader io.Reader, dbname string) TransAnnoScanner {
	scanner := io.NewGenericsTSVScanner[TransAnno](reader)
	return TransAnnoScanner{GenericsTSVScanner: scanner, DBname: dbname}
}

func (this TransAnnoScanner) Row() (TransAnnoRow, error) {
	row, err := this.StringMapRow()
	if err != nil {
		return TransAnnoRow{}, err
	}
	return TransAnnoRow{
		Snv:    fmt.Sprintf("%s\t%s\t%s\t%s\t%s", row["Chr"], row["Start"], row["End"], row["Ref"], row["Alt"]),
		Gene:   row[this.DBname+".Gene"],
		GeneID: row[this.DBname+".GeneID"],
		Event:  row[this.DBname+".Event"],
		Region: row[this.DBname+".Region"],
		Detail: row[this.DBname+".Detail"],
	}, nil
}

func AggsTransAnno(transAnnoOutput, dbname, outfile string) error {
	reader, err := io.NewIoReader(transAnnoOutput)
	if err != nil {
		return err
	}
	defer reader.Close()
	scanner := NewTransAnnoScanner(reader, dbname)
	geneAnnos := make(map[string]GeneAnno)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return err
		}
		geneAnno, ok := geneAnnos[row.ID()]
		if !ok {
			geneAnno = GeneAnno{
				Snv:    row.Snv,
				Gene:   row.Gene,
				GeneID: row.GeneID,
			}
		}
		geneAnno.AddInfo(row.Event, row.Region, row.Detail)
		geneAnnos[row.ID()] = geneAnno
	}
	geneAnnoTexts := make([]string, len(geneAnnos))
	i := 0
	for _, snvGeneAnno := range geneAnnos {
		geneAnnoTexts[i] = snvGeneAnno.Text()
		i++
	}
	sort.Strings(geneAnnoTexts)
	writer, err := os.Create(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprintf(writer, "%s\n", strings.Join(scanner.FieldNames, "\t"))
	for _, text := range geneAnnoTexts {
		fmt.Fprintf(writer, "%s\n", text)
	}
	return nil
}
