package core

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
)

type Position struct {
	ExonStart  int
	ExonEnd    int
	CdsStart   int
	CdsEnd     int
	ExonStarts []int
	ExonEnds   []int
}

type Region struct {
	Start     int
	End       int
	Typo      string
	ExonOrder int
}

type Regions []Region

type Refgene struct {
	Chrom      string
	Strand     byte
	Gene       string
	EntrezId   int
	Transcript string
	Position   Position
	Regions    Regions
	Streams    Regions
	Tag        string
	Mrna       Sequence
	Cdna       Sequence
	Protein    Sequence
}

type Refgenes []Refgene

type RefgeneDict map[string]Refgene

func (regions Regions) GetPrev(currentIndex int, strand byte) (Region, bool) {
	var index int
	if strand == '+' {
		index = currentIndex - 1
	} else {
		index = currentIndex + 1
	}
	if index >= 0 && index < len(regions) {
		return regions[index], true
	}
	return Region{}, false
}

func (regions Regions) GetNext(currentIndex int, strand byte) (Region, bool) {
	var index int
	if strand == '+' {
		index = currentIndex + 1
	} else {
		index = currentIndex - 1
	}
	if index >= 0 && index < len(regions) {
		return regions[index], true
	}
	return Region{}, false
}

func (regions Regions) Len() int {
	return len(regions)
}

func (regions Regions) Less(i, j int) bool {
	return regions[i].Start < regions[j].Start
}

func (regions Regions) Swap(i, j int) {
	regions[i], regions[j] = regions[j], regions[i]
}

func (position *Position) Read(refgeneLineField []string) {
	var err error
	if position.ExonStart, err = strconv.Atoi(refgeneLineField[4]); err != nil {
		panic(err)
	} else {
		position.ExonStart += 1
	}
	if position.ExonEnd, err = strconv.Atoi(refgeneLineField[5]); err != nil {
		panic(err)
	}
	if position.CdsStart, err = strconv.Atoi(refgeneLineField[6]); err != nil {
		panic(err)
	} else {
		position.CdsStart += 1
	}
	if position.CdsEnd, err = strconv.Atoi(refgeneLineField[7]); err != nil {
		panic(err)
	}
	for _, strPos := range strings.Split(strings.Trim(refgeneLineField[9], ","), ",") {
		if pos, err := strconv.Atoi(strPos); err == nil {
			position.ExonStarts = append(position.ExonStarts, pos+1)
		}
	}
	for _, strPos := range strings.Split(strings.Trim(refgeneLineField[10], ","), ",") {
		if pos, err := strconv.Atoi(strPos); err == nil {
			position.ExonEnds = append(position.ExonEnds, pos)
		}
	}
}

func (refgene Refgene) GetSn() string {
	return fmt.Sprintf("%s|%s:%d:%d", refgene.Transcript, refgene.Chrom, refgene.Position.ExonStart, refgene.Position.ExonEnd)
}

func (refgene Refgene) IsCmpl() bool {
	return refgene.Tag == "cmpl"
}

func (refgene *Refgene) SetUpDownStream(upDownStreamLen int) {
	stream1 := Region{
		Start: refgene.Position.ExonStart - upDownStreamLen,
		End:   refgene.Position.ExonStart - 1,
	}
	stream2 := Region{
		Start: refgene.Position.ExonEnd + 1,
		End:   refgene.Position.ExonEnd + upDownStreamLen,
	}
	if refgene.Strand == '+' {
		stream1.Typo = "upstream"
		stream2.Typo = "downstream"
	} else {
		stream1.Typo = "upstream"
		stream2.Typo = "downstream"
	}
	refgene.Streams = Regions{stream1, stream2}
}

func (refgene *Refgene) SetSequence(sequence Sequence) {
	if !sequence.IsEmpty() {
		refgene.Mrna = sequence
		if refgene.Tag != "unk" {
			for _, region := range refgene.Regions {
				if region.Typo == "cds" {
					seq := refgene.Mrna.GetSeq(region.Start-refgene.Position.ExonStart, region.End-region.Start+1)
					refgene.Cdna.Join([]Sequence{seq})
				}
			}
			if refgene.Strand == '-' {
				refgene.Cdna.Reverse()
			}
		}
		if !refgene.Cdna.IsEmpty() {
			refgene.Protein = refgene.Cdna.Translate(refgene.Chrom == "MT")
			if refgene.Protein.IsCmpl() {
				refgene.Tag = "cmpl"
			} else {
				refgene.Tag = "incmpl"
				refgene.Protein.Clear()
			}
		}
	}
}

func (refgene *Refgene) Read(refgeneLine string, ncbiGene NcbiGene) {
	field := strings.Split(refgeneLine, "\t")
	refgene.Position.Read(field)
	refgene.Chrom = strings.Replace(field[2], "chr", "", 1)
	refgene.Transcript = field[1]
	refgene.Strand = field[3][0]
	refgene.Gene = field[12]
	refgene.Tag = field[13]
	refgene.EntrezId = ncbiGene.GetEntrezId(refgene.Gene)
	exonNum := len(refgene.Position.ExonStarts)
	for i := 0; i < exonNum; i++ {
		if i > 0 {
			region := Region{
				Start: refgene.Position.ExonEnds[i-1] + 1,
				End:   refgene.Position.ExonStarts[i] - 1,
				Typo:  "intron",
			}
			refgene.Regions = append(refgene.Regions, region)
		}
		var exonOrder int
		if refgene.Strand == '+' {
			exonOrder = i + 1
		} else {
			exonOrder = exonNum - 1
		}
		start, end := refgene.Position.ExonStarts[i], refgene.Position.ExonEnds[i]
		if refgene.Position.CdsStart > end || refgene.Position.CdsEnd < start ||
			refgene.Position.CdsStart <= start && end <= refgene.Position.CdsEnd {
			var typo string
			if refgene.Position.CdsStart > end {
				if refgene.Strand == '+' {
					typo = "utr5"
				} else {
					typo = "utr3"
				}
			} else if refgene.Position.CdsEnd < start {
				if refgene.Strand == '-' {
					typo = "utr5"
				} else {
					typo = "utr3"
				}
			} else {
				typo = "cds"
			}
			refgene.Regions = append(refgene.Regions, Region{
				Start:     start,
				End:       end,
				Typo:      typo,
				ExonOrder: exonOrder,
			})
		} else {
			utrTypo1, utrTypo2 := "", ""
			cdsStart, cdsEnd := start, end
			if start < refgene.Position.CdsStart && refgene.Position.CdsStart < end {
				if refgene.Strand == '+' {
					utrTypo1 = "utr5"
				} else {
					utrTypo2 = "utr3"
				}
				cdsStart = refgene.Position.CdsStart
			}
			if start < refgene.Position.CdsEnd && refgene.Position.CdsEnd < end {
				if refgene.Strand == '+' {
					utrTypo2 = "utr3"
				} else {
					utrTypo2 = "utr5"
				}
				cdsEnd = refgene.Position.CdsEnd
			}
			if utrTypo1 != "" {
				refgene.Regions = append(refgene.Regions, Region{
					Start:     start,
					End:       refgene.Position.CdsStart - 1,
					Typo:      utrTypo1,
					ExonOrder: exonOrder,
				})
			}
			refgene.Regions = append(refgene.Regions, Region{
				Start:     cdsStart,
				End:       cdsEnd,
				Typo:      "cds",
				ExonOrder: exonOrder,
			})
			if utrTypo2 != "" {
				refgene.Regions = append(refgene.Regions, Region{
					Start:     refgene.Position.CdsEnd + 1,
					End:       end,
					Typo:      utrTypo2,
					ExonOrder: exonOrder,
				})
			}
		}
	}
	sort.Sort(refgene.Regions)
}

func (refgeneDict RefgeneDict) Read(refgeneFile string, ncbiGene NcbiGene) {
	if fp, err := os.Open(refgeneFile); err == nil {
		defer fp.Close()
		reader := bufio.NewReader(fp)
		for {
			if line, err := reader.ReadString('\n'); err == nil {
				line = strings.TrimSpace(line)
				if len(line) == 0 || line[0] == '#' {
					continue
				}
				var refgene Refgene
				refgene.Read(line, ncbiGene)
				if refgene.Chrom == "M" || len(refgene.Chrom) > 2 {
					continue
				}
				refgeneDict[refgene.GetSn()] = refgene
			} else {
				if err == io.EOF {
					break
				} else {
					panic(err)
				}
			}
		}
	} else {
		panic(err)
	}
}

func (refgeneDict RefgeneDict) SetUpDownStream(upDownStreamLen int) {
	for sn, refgene := range refgeneDict {
		refgene.SetUpDownStream(upDownStreamLen)
		refgeneDict[sn] = refgene
	}
}

func (refgeneDict RefgeneDict) SetSequence(mrna Fasta) {
	for sn, refgene := range refgeneDict {
		if sequence, ok := mrna[sn]; ok {
			refgene.SetSequence(sequence)
			refgeneDict[sn] = refgene
		}
	}
}
