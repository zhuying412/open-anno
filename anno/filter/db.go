package db

import (
	"bytes"
	"fmt"
	"io"
	"open-anno/anno/variant"
	"open-anno/pkg"
	"os"
	"strconv"
	"strings"
)

// FilterVarIdx FilterBased Variant Index
type FilterVarIdx struct {
	Bin   string `json:"bin"`
	Start int64  `json:"start"`
	End   int64  `json:"end"`
}

// ReadFilterVarIdx 读取 FilterBased Variant Index 文件
func ReadFilterVarIdx(infile string) (map[string]FilterVarIdx, int, error) {
	idxMap := make(map[string]FilterVarIdx)
	var binSize int
	reader, err := pkg.NewIOReader(infile)
	if err != nil {
		return idxMap, binSize, err
	}
	scanner := pkg.NewIOScanner(reader)
	scanner.Scan()
	binSize, err = strconv.Atoi(strings.Split(scanner.Text(), "\t")[1])
	if err != nil {
		return idxMap, binSize, err
	}
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		idx := FilterVarIdx{Bin: fmt.Sprintf("%s\t%s", row[0], row[1])}
		idx.Start, err = strconv.ParseInt(row[2], 10, 0)
		if err != nil {
			return idxMap, binSize, err
		}
		idx.End, err = strconv.ParseInt(row[3], 10, 0)
		if err != nil {
			return idxMap, binSize, err
		}
		idxMap[idx.Bin] = idx
	}
	return idxMap, binSize, err
}

// FilterVar FilterBased Variant
type FilterVariant struct {
	variant.AnnoVariant
	Info map[string]string `json:"info"`
}

// FilterVarReader FilterBased Variant Reader
type FilterVarReader struct {
	*os.File
	Header []string
}

// ReadRows 按Idx读取部分Variant
func (this FilterVarReader) ReadRows(idx FilterVarIdx) ([]FilterVariant, error) {
	this.File.Seek(idx.Start, io.SeekStart)
	buffer := make([]byte, idx.End-idx.Start)
	_, err := this.Read(buffer)
	if err != nil {
		return []FilterVariant{}, err
	}
	byteLines := bytes.Split(bytes.TrimSpace(buffer), []byte{'\n'})
	variants := make([]FilterVariant, len(byteLines))
	for i, byteLine := range byteLines {
		line := string(byteLine)
		row := strings.Split(line, "\t")
		variant := FilterVariant{
			AnnoVariant: variant.AnnoVariant{
				Chrom: row[0],
				Ref:   row[3],
				Alt:   row[4],
			},
		}
		variant.Start, err = strconv.Atoi(row[1])
		if err != nil {
			return variants, err
		}
		variant.End, err = strconv.Atoi(row[2])
		if err != nil {
			return variants, err
		}
		info := make(map[string]string)
		for i, header := range this.Header {
			info[header] = row[5+i]
		}
		variant.Info = info
		variants[i] = variant
	}
	return variants, nil
}

// FilterVarReader 新建FilterBased Variant Reader
func NewFilterVarReader(infile string) (FilterVarReader, error) {
	reader, err := os.Open(infile)
	if err != nil {
		return FilterVarReader{}, err
	}
	scanner := pkg.NewIOScanner(reader)
	scanner.Scan()
	return FilterVarReader{File: reader, Header: strings.Split(strings.TrimLeft(scanner.Text(), "#"), "\t")}, nil
}
