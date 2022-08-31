package io

import (
	"bytes"
	"fmt"
	"strings"
)

type AnnoResult struct {
	ID   string `json:"id"`
	Text string `json:"text"`
}

type AnnoResultScanner struct {
	Scanner[AnnoResult]
	FieldNames []string
}

func (this AnnoResultScanner) Header() string {
	return strings.Join(this.FieldNames, "\t")
}

func (this AnnoResultScanner) FillDot() string {
	var buffer bytes.Buffer
	for i, _ := range this.FieldNames {
		if i == 0 {
			buffer.WriteString(".")
		} else {
			buffer.WriteString("\t.")
		}
	}
	return buffer.String()
}

func NewAnnoResultScanner(reader Reader) AnnoResultScanner {
	scanner := NewScanner[AnnoResult](reader)
	scanner.Scan()
	fieldNames := strings.Split(scanner.Text(), "\t")[5:]
	return AnnoResultScanner{Scanner: scanner, FieldNames: fieldNames}
}

func (this AnnoResultScanner) Row() AnnoResult {
	fields := strings.Split(this.Text(), "\t")
	return AnnoResult{
		ID:   strings.Join(fields[0:5], ":"),
		Text: strings.Join(fields[5:], "\t"),
	}
}

func MergeAnnoResult(outfile, annoInput string, annoOuputs ...string) error {
	writer, err := NewIoWriter(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprintf(writer, "Chr\tStart\tEnd\tRef\tAlt\t")
	scanners := make([]AnnoResultScanner, len(annoOuputs))
	for i, annoOutput := range annoOuputs {
		reader, err := NewIoReader(annoOutput)
		if err != nil {
			return err
		}
		defer reader.Close()
		scanner := NewAnnoResultScanner(reader)
		scanners[i] = scanner
		fmt.Fprintf(writer, strings.Join(scanner.FieldNames, "\t"))
	}
	fmt.Fprint(writer, "\tOtherInfo\n")
	results := make([]map[string]string, len(scanners))
	for _, scanner := range scanners {
		result := make(map[string]string)
		for scanner.Scan() {
			row := scanner.Row()
			result[row.ID] = row.Text
		}
		results = append(results, result)
	}
	reader, err := NewIoReader(annoInput)
	if err != nil {
		return err
	}
	defer reader.Close()
	scanner := NewAnnoResultScanner(reader)
	for scanner.Scan() {
		row := scanner.Row()
		fmt.Fprint(writer, row.ID)
		for i, result := range results {
			if text, ok := result[row.ID]; ok {
				fmt.Fprintf(writer, "\t%s", text)
			} else {
				fmt.Fprintf(writer, "\t%s", scanners[i].FillDot())
			}
		}
		fmt.Fprintf(writer, "\t%s\n", row.Text)
	}
	return nil
}
