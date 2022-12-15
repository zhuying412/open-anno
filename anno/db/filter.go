package db

import (
	"bytes"
	"fmt"
	"io"
	"log"
	"open-anno/anno"
	"open-anno/pkg"
	"os"
	"sort"
	"strconv"
	"strings"
)

// CurBin 计算Bin游标
func CurBin[T int | int64](chrom string, start T, size int) string {
	return fmt.Sprintf("%s\t%d", chrom, start-(start%T(size)))
}

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

// FilterVarReader FilterBased Variant Reader
type FilterVarReader struct {
	*os.File
	Header string
}

// ReadRows 按Idx读取部分Variant
func (this FilterVarReader) ReadRows(idx FilterVarIdx) ([]anno.Variant, error) {
	this.File.Seek(idx.Start, io.SeekStart)
	buffer := make([]byte, idx.End-idx.Start)
	_, err := this.Read(buffer)
	if err != nil {
		return []anno.Variant{}, err
	}
	byteLines := bytes.Split(bytes.TrimSpace(buffer), []byte{'\n'})
	variants := make([]anno.Variant, len(byteLines))
	for i, byteLine := range byteLines {
		line := string(byteLine)
		row := strings.Split(line, "\t")
		variant := anno.Variant{
			Chrom:     row[0],
			Ref:       row[3],
			Alt:       row[4],
			Otherinfo: line,
		}
		variant.Start, err = strconv.Atoi(row[1])
		if err != nil {
			return variants, err
		}
		variant.End, err = strconv.Atoi(row[2])
		if err != nil {
			return variants, err
		}
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
	return FilterVarReader{File: reader, Header: strings.TrimLeft(scanner.Text(), "#")}, nil
}

// ReadAnnoInput 读取AnnoInput，并按Bin游标分类保存
func ReadAnnoInput(annoInputFile string, binSize int) (map[string]anno.Variants, error) {
	variants := make(map[string]anno.Variants)
	reader, err := pkg.NewIOReader(annoInputFile)
	if err != nil {
		return variants, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		variant := anno.Variant{Chrom: row[0], Ref: row[3], Alt: row[4], Otherinfo: row[5]}
		variant.Start, err = strconv.Atoi(row[1])
		if err != nil {
			return variants, err
		}
		variant.End, err = strconv.Atoi(row[2])
		if err != nil {
			return variants, err
		}
		curbin := CurBin(variant.Chrom, variant.Start, binSize)
		if _, ok := variants[curbin]; !ok {
			variants[curbin] = make(anno.Variants, 0)
		}
		variants[curbin] = append(variants[curbin], variant)
	}
	return variants, err
}

// AnnoFilterBased 注释FilterBased
func AnnoFilterBased(annoInputFile, dbFile, dbIdxFile, outfile string) error {
	// 读取DB Idx File输入文件
	log.Printf("Read DBIdxFile: %s ...", dbIdxFile)
	idxMap, binSize, err := ReadFilterVarIdx(dbIdxFile)
	if err != nil {
		return err
	}
	// 读取变异输入文件
	log.Printf("Read AnnoInput: %s ...", annoInputFile)
	variantMap, err := ReadAnnoInput(annoInputFile, binSize)
	if err != nil {
		return err
	}
	// 定义输出
	writer, err := pkg.NewIOWriter(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	// 开始注释
	reader, err := NewFilterVarReader(dbFile)
	if err != nil {
		return err
	}
	defer reader.Close()
	fmt.Fprintf(writer, "%s\n", reader.Header)
	for curbin, variants := range variantMap {
		sort.Sort(variants)
		if idx, ok := idxMap[curbin]; ok {
			dbvars, err := reader.ReadRows(idx)
			if err != nil {
				return err
			}
			for i, j := 0, 0; i < len(variants) && j < len(dbvars); {
				if dbvars[j].Greater(variants[i]) {
					i++
				} else if dbvars[j].Less(variants[i]) {
					j++
				} else {
					fmt.Fprintf(writer, "%s\n", dbvars[j].Otherinfo)
					i++
				}
			}
		}
	}
	return nil
}
