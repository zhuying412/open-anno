package refgene

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"open-anno/pkg/seq"
	"strconv"
	"strings"

	"github.com/brentp/faidx"
)

// Transcript Refgene
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
	GeneID     string  `json:"gene_id"`
	CdsStat    string  `json:"cdsStartStat"`
	Regions    Regions `json:"regions"`
}

func (this Transcript) CdsCount() int {
	var count int
	for _, region := range this.Regions {
		if region.Type == RType_CDS {
			count++
		}
	}
	return count
}

func (this Transcript) CLen() int {
	var length int
	for _, region := range this.Regions {
		if region.Type == RType_CDS {
			length += region.End - region.Start + 1
		}
	}
	return length
}

func (this Transcript) CDNA() string {
	var cdna bytes.Buffer
	for _, region := range this.Regions {
		if region.Type == RType_CDS {
			cdna.WriteString(region.Sequence)
		}
	}
	return cdna.String()
}

func (this Transcript) DNA() string {
	var cdna bytes.Buffer
	for _, region := range this.Regions {
		cdna.WriteString(region.Sequence)
	}
	return cdna.String()
}

func (this Transcript) IsUnk() bool {
	return this.CdsEnd-this.CdsStart+1 == 0
}

func (this *Transcript) SetRegions(mrna *faidx.Faidx, symbolToId map[string]string, seqRequired bool) error {
	regions, err := NewRegions(*this)
	if err != nil {
		return err
	}
	for i, region := range regions {
		mrnaName := fmt.Sprintf("%s:%s:%s", this.Chrom, this.Gene, this.Name)
		regions[i].Sequence, err = seq.Fetch(mrna, mrnaName, region.Start-this.TxStart, region.End-this.TxStart+1)
		if err != nil && seqRequired {
			return err
		}
	}
	this.Regions = regions
	this.GeneID = symbolToId[this.Gene]
	return nil
}

// Transcripts
type Transcripts map[string]Transcript

func (this Transcripts) FilterChrom(chrom string, symbolToId map[string]string, mrna *faidx.Faidx, seqRequired bool) (Transcripts, error) {
	transcripts := make(Transcripts)
	for sn, trans := range this {
		if trans.Chrom == chrom {
			err := trans.SetRegions(mrna, symbolToId, seqRequired)
			if err != nil {
				return transcripts, err
			}
			transcripts[sn] = trans
		}
	}
	return transcripts, nil
}

type TransScanner struct {
	scanner *bufio.Scanner
}

func NewTransScanner(reader io.Reader) TransScanner {
	scanner := bufio.NewScanner(reader)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)
	return TransScanner{scanner: scanner}
}

func (this *TransScanner) Scan() bool {
	return this.scanner.Scan()
}

func (this TransScanner) Text() string {
	return this.scanner.Text()
}

func (this TransScanner) Row() (Transcript, error) {
	transcript := Transcript{}
	fields := strings.Split(this.Text(), "\t")
	var name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, gene string
	var exonStarts, exonEnds []string
	switch len(fields) {
	case 16:
		name = fields[1]
		chrom = pkg.FormatChrom(fields[2])
		strand = fields[3]
		txStart = fields[4]
		txEnd = fields[5]
		cdsStart = fields[6]
		cdsEnd = fields[7]
		exonCount = fields[8]
		exonStarts = strings.Split(strings.Trim(fields[9], ","), ",")
		exonEnds = strings.Split(strings.Trim(fields[10], ","), ",")
		gene = fields[12]
	case 12:
		name = fields[0]
		chrom = pkg.FormatChrom(fields[1])
		strand = fields[2]
		txStart = fields[3]
		txEnd = fields[4]
		cdsStart = fields[5]
		cdsEnd = fields[6]
		exonCount = fields[7]
		exonStarts = strings.Split(strings.Trim(fields[8], ","), ",")
		exonEnds = strings.Split(strings.Trim(fields[9], ","), ",")
		gene = fields[0]
	default:
		return transcript, errors.New("unkown refGene format")
	}
	var err error
	transcript.Name = name
	transcript.Chrom = chrom
	transcript.Strand = strand
	transcript.Gene = gene
	transcript.TxStart, err = strconv.Atoi(txStart)
	if err != nil {
		return transcript, err
	}
	transcript.TxStart++
	transcript.TxEnd, err = strconv.Atoi(txEnd)
	if err != nil {
		return transcript, err
	}
	transcript.CdsStart, err = strconv.Atoi(cdsStart)
	if err != nil {
		return transcript, err
	}
	transcript.CdsStart++
	transcript.CdsEnd, err = strconv.Atoi(cdsEnd)
	if err != nil {
		return transcript, err
	}
	transcript.ExonCount, err = strconv.Atoi(exonCount)
	if err != nil {
		return transcript, err
	}
	transcript.ExonStarts = make([]int, transcript.ExonCount)
	transcript.ExonEnds = make([]int, transcript.ExonCount)
	for i := 0; i < transcript.ExonCount; i++ {
		transcript.ExonStarts[i], err = strconv.Atoi(exonStarts[i])
		if err != nil {
			return transcript, err
		}
		transcript.ExonStarts[i]++
		transcript.ExonEnds[i], err = strconv.Atoi(exonEnds[i])
		if err != nil {
			return transcript, err
		}
	}
	return transcript, err
}

func ReadRefgene(infile string) (Transcripts, error) {
	transcripts := make(Transcripts)
	fi, err := io.NewIoReader(infile)
	if err != nil {
		return transcripts, err
	}
	defer fi.Close()
	scanner := NewTransScanner(fi)
	for scanner.Scan() {
		transcript, err := scanner.Row()
		if err != nil {
			return transcripts, err
		}
		if _, ok := seq.GENOME[transcript.Chrom]; ok {
			transcripts[transcript.Name] = transcript
		}
	}
	return transcripts, err
}
