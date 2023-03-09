package db

import (
	"bytes"
	"fmt"
	"io"
	"open-anno/anno"
	"open-anno/pkg"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/brentp/vcfgo"
)

// FilterVarIdx FilterBased Variant Index
type FilterVarIdx struct {
	Bin   string `json:"bin"`
	Start int64  `json:"start"`
	End   int64  `json:"end"`
}

// ReadFilterVarIdx 读取 FilterBased Variant Index 文件
func ReadFilterVarIdx(infile string) (map[string]FilterVarIdx, int, error) {
	idxs := make(map[string]FilterVarIdx)
	var binSize int
	reader, err := pkg.NewIOReader(infile)
	if err != nil {
		return idxs, binSize, err
	}
	scanner := pkg.NewIOScanner(reader)
	scanner.Scan()
	binSize, err = strconv.Atoi(strings.Split(scanner.Text(), "\t")[1])
	if err != nil {
		return idxs, binSize, err
	}
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		idx := FilterVarIdx{Bin: fmt.Sprintf("%s\t%s", row[0], row[1])}
		idx.Start, err = strconv.ParseInt(row[2], 10, 0)
		if err != nil {
			return idxs, binSize, err
		}
		idx.End, err = strconv.ParseInt(row[3], 10, 0)
		if err != nil {
			return idxs, binSize, err
		}
		idxs[idx.Bin] = idx
	}
	return idxs, binSize, err
}

// FilterVar FilterBased Variant
type FilterVariant struct {
	anno.AnnoVariant
	Info []anno.AnnoInfo `json:"info"`
}

// FilterVarReader FilterBased Variant Reader
type FilterVarReader struct {
	os.File
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
	filterVariants := make([]FilterVariant, len(byteLines))
	for i, byteLine := range byteLines {
		line := string(byteLine)
		row := strings.Split(line, "\t")
		filterVariant := FilterVariant{
			AnnoVariant: anno.AnnoVariant{
				Chrom: row[0],
				Ref:   row[3],
				Alt:   row[4],
			},
		}
		filterVariant.Start, err = strconv.Atoi(row[1])
		if err != nil {
			return filterVariants, err
		}
		filterVariant.End, err = strconv.Atoi(row[2])
		if err != nil {
			return filterVariants, err
		}
		infos := make([]anno.AnnoInfo, len(this.Header)-5)
		for i := 5; i < len(this.Header); i++ {
			infos[i-5] = anno.AnnoInfo{Key: this.Header[i], Value: row[i]}
		}
		filterVariant.Info = infos
		filterVariants[i] = filterVariant
	}
	return filterVariants, nil
}

func (this FilterVarReader) VCFHeaderInfo() map[string]*vcfgo.Info {
	info := make(map[string]*vcfgo.Info)
	for _, header := range this.Header {
		info[header] = &vcfgo.Info{
			Id:          header,
			Description: "Annotation of " + header,
			Number:      ".",
			Type:        "UNKNOWN",
		}
	}
	return info
}

// FilterVarReader 新建FilterBased Variant Reader
func NewFilterVarReader(infile string) (*FilterVarReader, error) {
	reader, err := os.Open(infile)
	if err != nil {
		return &FilterVarReader{}, err
	}
	scanner := pkg.NewIOScanner(reader)
	scanner.Scan()
	return &FilterVarReader{File: *reader, Header: strings.Split(strings.TrimLeft(scanner.Text(), "#"), "\t")}, nil
}

// AnnoFilterBased 注释SNV FilterBased
func AnnoFilterBased(variants anno.Variants, filterVarReader *FilterVarReader, filterVarIdxs map[string]FilterVarIdx, binSize int) (anno.AnnoInfos, error) {
	annoInfos := make(anno.AnnoInfos)
	// 读取DB Idx File输入文件
	// log.Printf("Read DBIdxFile: %s ...", dbIdxFile)
	// idxMap, binSize, err := ReadFilterVarIdx(dbIdxFile)
	// if err != nil {
	// 	return annoInfos, err
	// }
	for curbin, snvs := range variants.AggregateByBin(binSize) {
		sort.Sort(snvs)
		if idx, ok := filterVarIdxs[curbin]; ok {
			dbvars, err := filterVarReader.ReadRows(idx)
			if err != nil {
				return annoInfos, err
			}
			for i, j := 0, 0; i < len(snvs) && j < len(dbvars); {
				annoVariant := snvs[i].AnnoVariant()
				if dbvars[j].Greater(annoVariant) {
					i++
				} else if dbvars[j].Less(annoVariant) {
					j++
				} else {
					annoInfos[annoVariant.PK()] = dbvars[j].Info
					i++
				}
			}
		}
	}
	return annoInfos, nil
}
