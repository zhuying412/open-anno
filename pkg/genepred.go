package pkg

import (
	"errors"
	"fmt"
	"strconv"
	"strings"
)

// GenePred GenePred文件解析结果（如refGene等）
type GenePred struct {
	Name       string `json:"name"`
	Chrom      string `json:"chrom"`
	Strand     string `json:"strand"`
	TxStart    int    `json:"txStart"`
	TxEnd      int    `json:"txEnd"`
	CdsStart   int    `json:"cdsStart"`
	CdsEnd     int    `json:"cdsEnd"`
	ExonCount  int    `json:"exonCount"`
	ExonStarts []int  `json:"exonStarts"`
	ExonEnds   []int  `json:"exonEnds"`
	Gene       string `json:"gene"`
	CdsStat    string `json:"cdsStartStat"`
}

// PK GenePred主键名称
func (this GenePred) PK() string {
	return fmt.Sprintf("%s:%s:%s", this.Chrom, this.Gene, this.Name)
}

// IsUnk 是否为ncRNA
func (this GenePred) IsUnk() bool {
	return this.CdsEnd-this.CdsStart+1 == 0
}

// NewGenePred 从GenePred文件中读取一行并解析为GenePred对象
func NewGenePred(row []string) (GenePred, error) {
	var gpe GenePred
	var name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, gene string
	var exonStarts, exonEnds []string
	switch len(row) {
	case 16:
		name = row[1]
		chrom = row[2]
		strand = row[3]
		txStart = row[4]
		txEnd = row[5]
		cdsStart = row[6]
		cdsEnd = row[7]
		exonCount = row[8]
		exonStarts = strings.Split(strings.Trim(row[9], ","), ",")
		exonEnds = strings.Split(strings.Trim(row[10], ","), ",")
		gene = row[12]
	case 12:
		name = row[0]
		chrom = row[1]
		strand = row[2]
		txStart = row[3]
		txEnd = row[4]
		cdsStart = row[5]
		cdsEnd = row[6]
		exonCount = row[7]
		exonStarts = strings.Split(strings.Trim(row[8], ","), ",")
		exonEnds = strings.Split(strings.Trim(row[9], ","), ",")
		gene = row[0]
	default:
		return gpe, errors.New("unknown refGene format")
	}
	var err error
	gpe.Name = name
	gpe.Chrom = chrom
	gpe.Strand = strand
	gpe.Gene = gene
	gpe.TxStart, err = strconv.Atoi(txStart)
	if err != nil {
		return gpe, err
	}
	gpe.TxStart++
	gpe.TxEnd, err = strconv.Atoi(txEnd)
	if err != nil {
		return gpe, err
	}
	gpe.CdsStart, err = strconv.Atoi(cdsStart)
	if err != nil {
		return gpe, err
	}
	gpe.CdsStart++
	gpe.CdsEnd, err = strconv.Atoi(cdsEnd)
	if err != nil {
		return gpe, err
	}
	gpe.ExonCount, err = strconv.Atoi(exonCount)
	if err != nil {
		return gpe, err
	}
	gpe.ExonStarts = make([]int, gpe.ExonCount)
	gpe.ExonEnds = make([]int, gpe.ExonCount)
	for i := 0; i < gpe.ExonCount; i++ {
		gpe.ExonStarts[i], err = strconv.Atoi(exonStarts[i])
		if err != nil {
			return gpe, err
		}
		gpe.ExonStarts[i]++
		gpe.ExonEnds[i], err = strconv.Atoi(exonEnds[i])
		if err != nil {
			return gpe, err
		}
	}
	return gpe, err
}
