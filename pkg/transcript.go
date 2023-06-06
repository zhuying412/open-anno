package pkg

import (
	"errors"
	"fmt"
	"strconv"
	"strings"

	"github.com/brentp/faidx"
)

var GeneSymbolToID map[string]map[string]string

func InitGeneSymbolToID(infile string) error {
	gene := make(map[string]map[string]string)
	reader, err := NewIOReader(infile)
	if err != nil {
		return err
	}
	defer reader.Close()
	scanner := NewCSVScanner(reader)
	scanner.Scan()
	for scanner.Scan() {
		row := scanner.Row()
		chrom := row["Chrom"]
		entrezId := row["EntrezId"]
		symbol := row["Symbol"]
		if _, ok := gene[chrom]; !ok {
			gene[chrom] = make(map[string]string)
		}
		gene[chrom][symbol] = entrezId

	}
	GeneSymbolToID = gene
	return nil
}

// Transcript 转录本，继承自GenePred，加入GeneID和Regions信息
type Transcript struct {
	Name       string  `json:"name"`
	Chrom      string  `json:"chrom"`
	Strand     string  `json:"strand"`
	TxStart    int     `json:"txStart"`
	TxEnd      int     `json:"txEnd"`
	CdsStart   int     `json:"cdsStart"`
	CdsEnd     int     `json:"cdsEnd"`
	ExonCount  int     `json:"exonCount"`
	ExonStarts []int   `json:"exonStarts"`
	ExonEnds   []int   `json:"exonEnds"`
	Gene       string  `json:"gene"`
	CdsStat    string  `json:"cdsStartStat"`
	GeneID     string  `json:"gene_id"`
	Regions    Regions `json:"regions"`
}

// PK GenePred主键名称
func (this Transcript) PK() string {
	return fmt.Sprintf("%s:%s:%s", this.Chrom, this.Gene, this.Name)
}

// IsUnk 是否为ncRNA
func (this Transcript) IsUnk() bool {
	return this.CdsEnd-this.CdsStart+1 == 0
}

// HasUTR3 存在UTR3区域
func (this Transcript) HasUTR3() bool {
	if this.Strand == "+" {
		return this.CdsEnd < this.TxEnd
	}
	return this.CdsStart > this.TxStart
}

// HasUTR5 存在UTR5区域
func (this Transcript) HasUTR5() bool {
	if this.Strand == "+" {
		return this.CdsStart > this.TxStart
	}
	return this.CdsEnd < this.TxEnd
}

// CdsCount Regions中CDS元件数量
func (this Transcript) CdsCount() int {
	return this.Regions.CdsCount()
}

// CLen Regions中CDS区域总长度
func (this Transcript) CLen() int {
	return this.Regions.CLen()
}

// ULen 转录本exon区UTR的长度
func (this Transcript) ULen() (int, int) {
	var utrLen1, utrLen2 int
	for _, region := range this.Regions {
		if region.Type == RType_UTR {
			if region.End < this.CdsStart {
				utrLen1 += region.Len()
			} else {
				utrLen2 += region.Len()
			}
		}
	}
	return utrLen1, utrLen2
}

// CDNA Regions中CDS的Sequence拼接结果即CDNA
func (this Transcript) CDNA() string {
	return this.Regions.CDNA()
}

// CDNA Regions中所有Sequence拼接结果即DNA
func (this Transcript) DNA() string {
	return this.Regions.DNA()
}

// Region 指定pos所在的region信息
func (this Transcript) Region(pos int) (Region, int, int) {
	var region Region
	var cLen, uLen int
	for _, region := range this.Regions {
		if region.Start <= pos && pos <= region.End {
			return region, cLen, uLen
		}
		if region.Type == RType_CDS {
			cLen += region.Len()
			uLen = 0
		}
		if region.Type == RType_UTR {
			uLen += region.Len()
		}
	}
	return region, cLen, uLen
}

// SetGeneID 根据geneSymbolToID的Map信息设置转录本的GeneID
func (this *Transcript) SetGeneID() {
	this.GeneID = "."
	if entrezId, ok := GeneSymbolToID[this.Chrom][this.Gene]; ok {
		this.GeneID = entrezId
	}
}

// SetRegions 设置初始化Regions信息
func (this *Transcript) SetRegions() {
	this.Regions = NewRegions(*this)
}

// SetRegions 设置初始化Regions信息, 并设置每个Regions的Sequence
func (this *Transcript) SetRegionsWithSeq(genome *faidx.Faidx) error {
	var err error
	seq, err := genome.Get(this.Chrom, this.TxStart-1, this.TxEnd)
	if err != nil {
		return err
	}
	this.Regions = NewRegionsWithSeq(*this, seq)
	return nil
}

// NewGenePred 从GenePred文件中读取一行并解析为GenePred对象
func NewTranscript(line string) (Transcript, error) {
	row := strings.Split(line, "\t")
	var trans Transcript
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
		return trans, errors.New("unknown refGene format")
	}
	var err error
	trans.Name = name
	trans.Chrom = chrom
	trans.Strand = strand
	trans.Gene = gene
	trans.TxStart, err = strconv.Atoi(txStart)
	if err != nil {
		return trans, err
	}
	trans.TxStart++
	trans.TxEnd, err = strconv.Atoi(txEnd)
	if err != nil {
		return trans, err
	}
	trans.CdsStart, err = strconv.Atoi(cdsStart)
	if err != nil {
		return trans, err
	}
	trans.CdsStart++
	trans.CdsEnd, err = strconv.Atoi(cdsEnd)
	if err != nil {
		return trans, err
	}
	trans.ExonCount, err = strconv.Atoi(exonCount)
	if err != nil {
		return trans, err
	}
	trans.ExonStarts = make([]int, trans.ExonCount)
	trans.ExonEnds = make([]int, trans.ExonCount)
	for i := 0; i < trans.ExonCount; i++ {
		trans.ExonStarts[i], err = strconv.Atoi(exonStarts[i])
		if err != nil {
			return trans, err
		}
		trans.ExonStarts[i]++
		trans.ExonEnds[i], err = strconv.Atoi(exonEnds[i])
		if err != nil {
			return trans, err
		}
	}
	return trans, err
}
