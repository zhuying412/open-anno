package gene

import (
	"bufio"
	"bytes"
	"io"
	"open-anno/pkg"
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

func (this *Transcript) SetRegions(fai *faidx.Faidx) error {
	var err error
	this.Regions, err = NewRegions(*this, fai)
	return err
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

// Transcripts
type Transcripts map[string]Transcript

func ReadRefgene(refgeneFile string) (Transcripts, error) {
	transcripts := make(Transcripts)
	fi, err := os.Open(refgeneFile)
	if err != nil {
		return transcripts, err
	}
	defer fi.Close()
	scanner := bufio.NewScanner(fi)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		chrom := strings.Replace(fields[2], "chr", "", -1)
		if _, ok := GENOME[chrom]; !ok {
			continue
		}
		transcript := Transcript{
			Name:         fields[1],
			Chrom:        chrom,
			Strand:       fields[3],
			Gene:         fields[12],
			CdsStartStat: fields[13],
			CdsEndStat:   fields[14],
		}
		transcript.Bin, err = strconv.Atoi(fields[0])
		if err != nil {
			return transcripts, err
		}
		transcript.TxStart, err = strconv.Atoi(fields[4])
		if err != nil {
			return transcripts, err
		}
		transcript.TxStart++
		transcript.TxEnd, err = strconv.Atoi(fields[5])
		if err != nil {
			return transcripts, err
		}
		transcript.CdsStart, err = strconv.Atoi(fields[6])
		if err != nil {
			return transcripts, err
		}
		transcript.CdsStart++
		transcript.CdsEnd, err = strconv.Atoi(fields[7])
		if err != nil {
			return transcripts, err
		}
		transcript.ExonCount, err = strconv.Atoi(fields[8])
		if err != nil {
			return transcripts, err
		}
		transcript.Score, err = strconv.Atoi(fields[11])
		if err != nil {
			return transcripts, err
		}
		for _, tpos := range strings.Split(strings.Trim(fields[9], ","), ",") {
			npos, err := strconv.Atoi(tpos)
			if err != nil {
				return transcripts, err
			}
			transcript.ExonStarts = append(transcript.ExonStarts, npos+1)
		}
		for _, tpos := range strings.Split(strings.Trim(fields[10], ","), ",") {
			npos, err := strconv.Atoi(tpos)
			if err != nil {
				return transcripts, err
			}
			transcript.ExonEnds = append(transcript.ExonEnds, npos)
		}
		for _, tpos := range strings.Split(strings.Trim(fields[15], ","), ",") {
			npos, err := strconv.Atoi(tpos)
			if err != nil {
				return transcripts, err
			}
			transcript.ExonFrames = append(transcript.ExonFrames, npos)
		}
		transcripts[transcript.Name] = transcript
	}
	return transcripts, err
}

func ReadTransDB(dbFile string) (Transcripts, error) {
	transcripts := make(Transcripts)
	fi, err := os.Open(dbFile)
	if err != nil {
		return transcripts, err
	}
	defer fi.Close()
	reader := bufio.NewReader(fi)
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			return transcripts, err
		}
		trans := pkg.FromJSON[Transcript](line)
		transcripts[trans.Name] = trans
	}
	return transcripts, err
}
