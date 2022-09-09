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
	for i := range this.FieldNames {
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

type AnnoResultScanners []AnnoResultScanner

func (this AnnoResultScanners) Header() string {
	fieldNames := make([]string, 0)
	for _, scanner := range this {
		fieldNames = append(fieldNames, scanner.FieldNames...)
	}
	return strings.Join(fieldNames, "\t")
}

func (this AnnoResultScanners) Results() []map[string]string {
	results := make([]map[string]string, len(this))
	for i, scanner := range this {
		result := make(map[string]string)
		for scanner.Scan() {
			row := scanner.Row()
			result[row.ID] = row.Text
		}
		results[i] = result
	}
	return results
}

func (this AnnoResultScanner) Row() AnnoResult {
	fields := strings.Split(this.Text(), "\t")
	return AnnoResult{
		ID:   strings.Join(fields[0:5], ":"),
		Text: strings.Join(fields[5:], "\t"),
	}
}

func MergeAnnoResult(outfile, annoInput, annoGBOutput string, annoOuputs ...string) error {
	writer, err := NewIoWriter(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	scanners := make(AnnoResultScanners, len(annoOuputs)+1)
	for i, annoOutput := range append(annoOuputs, annoInput) {
		reader, err := NewIoReader(annoOutput)
		if err != nil {
			return err
		}
		defer reader.Close()
		scanner := NewAnnoResultScanner(reader)
		scanners[i] = scanner
	}
	results := scanners.Results()
	reader, err := NewIoReader(annoGBOutput)
	if err != nil {
		return err
	}
	defer reader.Close()
	scanner := NewAnnoResultScanner(reader)
	fmt.Fprintf(writer, "Chr\tStart\tEnd\tRef\tAlt\t%s\t%s\n", scanner.Header(), scanners.Header())
	for scanner.Scan() {
		row := scanner.Row()
		fmt.Fprint(writer, scanner.Text())
		for i, result := range results {
			if text, ok := result[row.ID]; ok {
				fmt.Fprintf(writer, "\t%s", text)
			} else {
				fmt.Fprintf(writer, "\t%s", scanners[i].FillDot())
			}
		}
		fmt.Fprint(writer, "\n")
	}
	return nil
}
