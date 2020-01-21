package core

import (
	"bufio"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
)

type Refidx struct {
	Chrom       string
	Start       int
	End         int
	Transcripts []string
}

type Refidxs []Refidx

type RefidxDict map[string]Refidxs

func (refidx Refidx) GetDigitalPosition() (int, int) {
	start := ChromOrderDict[refidx.Chrom]*1e9 + refidx.Start
	end := ChromOrderDict[refidx.Chrom]*1e9 + refidx.End
	return start, end
}

func (refidx Refidx) GetRefgenes(refgeneDict RefgeneDict) Refgenes {
	var refgenes Refgenes
	for _, transcript := range refidx.Transcripts {
		if refgene, ok := refgeneDict[transcript]; ok {
			refgenes = append(refgenes, refgene)
		}
	}
	return refgenes
}

func (refidxs Refidxs) Len() int {
	return len(refidxs)
}

func (refidxs Refidxs) Less(i, j int) bool {
	starti, endi := refidxs[i].GetDigitalPosition()
	startj, endj := refidxs[j].GetDigitalPosition()
	if starti == startj {
		return endi < endj
	} else {
		return starti < startj
	}
}

func (refidxs Refidxs) Swap(i, j int) {
	refidxs[i], refidxs[j] = refidxs[j], refidxs[i]
}

func (refidxs *Refidxs) Read(refidxFile string) {
	if fp, err := os.Open(refidxFile); err == nil {
		defer fp.Close()
		reader := bufio.NewReader(fp)
		for {
			if line, err := reader.ReadString('\n'); err == nil {
				line = strings.TrimSpace(line)
				field := strings.Split(line, "\t")
				refidx := Refidx{
					Chrom:       field[0],
					Transcripts: strings.Split(field[3], ","),
				}
				if refidx.Start, err = strconv.Atoi(field[1]); err != nil {
					panic(err)
				}
				if refidx.End, err = strconv.Atoi(field[2]); err != nil {
					panic(err)
				}
				*refidxs = append(*refidxs, refidx)
			} else {
				if err == io.EOF {
					break
				} else {
					panic(err)
				}
			}
		}
		sort.Sort(refidxs)
	} else {
		panic(err)
	}
}
