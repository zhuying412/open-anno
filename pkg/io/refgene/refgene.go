package refgene

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/seq"
	"os"
	"strconv"
	"strings"

	"github.com/brentp/faidx"
)

const (
	CdsStat_N = "none"
	CdsStat_C = "cmpl"
	CdsStat_I = "incompl"
	CdsStat_U = "unk"
)

// Transcript Refgene
type Transcript struct {
	Bin          int     `json:"bin"`
	Name         string  `json:"name"`
	Chrom        string  `json:"chrom"`
	Strand       string  `json:"strand"`
	TxStart      int     `json:"txStart"`
	TxEnd        int     `json:"txEnd"`
	CdsStart     int     `json:"cdsStart"`
	CdsEnd       int     `json:"cdsEnd"`
	ExonCount    int     `json:"exonCount"`
	ExonStarts   []int   `json:"exonStarts"`
	ExonEnds     []int   `json:"exonEnds"`
	Score        int     `json:"score"`
	Gene         string  `json:"gene"`
	GeneID       string  `json:"gene_id"`
	CdsStartStat string  `json:"cdsStartStat"`
	CdsEndStat   string  `json:"cdsEndStat"`
	ExonFrames   []int   `json:"exonFrames"`
	Regions      Regions `json:"regions"`
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

func (this Transcript) IsCmpl() bool {
	return this.CdsStartStat == CdsStat_C && this.CdsEndStat == CdsStat_C
}

func (this Transcript) IsUnk() bool {
	return this.CdsStartStat == CdsStat_U && this.CdsEndStat == CdsStat_U
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
			log.Fatal(err)
		}
	}
	this.Regions = regions
	this.GeneID = symbolToId[this.Gene]
	return nil
}

// Transcripts
type Transcripts map[string]Transcript

func (this Transcripts) FilterChrom(chrom string, symbolToId map[string]string, mrna *faidx.Faidx, seqRequired bool) Transcripts {
	transcripts := make(Transcripts)
	for sn, trans := range this {
		if trans.Chrom == chrom {
			err := trans.SetRegions(mrna, symbolToId, seqRequired)
			if err != nil {
				log.Fatal(err)
			}
			transcripts[sn] = trans
		}
	}
	return transcripts
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
	fields := strings.Split(this.Text(), "\t")
	chrom := pkg.FormatChrom(fields[2])
	transcript := Transcript{
		Name:         fields[1],
		Chrom:        chrom,
		Strand:       fields[3],
		Gene:         fields[12],
		CdsStartStat: fields[13],
		CdsEndStat:   fields[14],
	}
	var err error
	transcript.Bin, err = strconv.Atoi(fields[0])
	if err != nil {
		return transcript, err
	}
	transcript.TxStart, err = strconv.Atoi(fields[4])
	if err != nil {
		return transcript, err
	}
	transcript.TxStart++
	transcript.TxEnd, err = strconv.Atoi(fields[5])
	if err != nil {
		return transcript, err
	}
	transcript.CdsStart, err = strconv.Atoi(fields[6])
	if err != nil {
		return transcript, err
	}
	transcript.CdsStart++
	transcript.CdsEnd, err = strconv.Atoi(fields[7])
	if err != nil {
		return transcript, err
	}
	transcript.ExonCount, err = strconv.Atoi(fields[8])
	if err != nil {
		return transcript, err
	}
	transcript.Score, err = strconv.Atoi(fields[11])
	if err != nil {
		return transcript, err
	}
	exonStarts := strings.Split(strings.Trim(fields[9], ","), ",")
	exonEnds := strings.Split(strings.Trim(fields[10], ","), ",")
	exonFrames := strings.Split(strings.Trim(fields[15], ","), ",")
	transcript.ExonStarts = make([]int, transcript.ExonCount)
	transcript.ExonEnds = make([]int, transcript.ExonCount)
	transcript.ExonFrames = make([]int, transcript.ExonCount)
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
		transcript.ExonFrames[i], err = strconv.Atoi(exonFrames[i])
		if err != nil {
			return transcript, err
		}
	}
	return transcript, err
}

func ReadRefgene(infile string) (Transcripts, error) {
	transcripts := make(Transcripts)
	fi, err := os.Open(infile)
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
		if _, ok := seq.GENOME[transcript.Chrom]; !ok {
			continue
		}
	}
	return transcripts, err
}
