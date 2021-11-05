package db

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

type Refgene struct {
	Chrom      string
	Start      int
	End        int
	UpStream   int
	DownStream int
	Transcript string
}

func NewRefgene(refgeneLine string) Refgene {
	refgene := Refgene{}
	field := strings.Split(refgeneLine, "\t")
	refgene.Transcript = field[1]
	refgene.Chrom = strings.Replace(strings.Split(field[2], " ")[0], "chr", "", -1)
	if start, err := strconv.Atoi(field[4]); err == nil {
		refgene.Start = start + 1
	} else {
		log.Panic(err)
	}
	if end, err := strconv.Atoi(field[5]); err == nil {
		refgene.End = end
	} else {
		log.Panic(err)
	}
	refgene.UpStream = refgene.Start - UpDownStreamLen
	refgene.DownStream = refgene.End + UpDownStreamLen
	return refgene
}

func (r Refgene) SN() string {
	return fmt.Sprintf("%s|%s:%d:%d", r.Transcript, r.Chrom, r.Start, r.End)
}

type Refgenes []Refgene

func (refgenes Refgenes) Len() int {
	return len(refgenes)
}

func (r Refgenes) Less(i, j int) bool {
	orderChromi, _ := ChromArray.GetByName(r[i].Chrom)
	orderChromj, _ := ChromArray.GetByName(r[j].Chrom)
	if orderChromi != orderChromj {
		return orderChromi < orderChromj
	}
	starti, endi := r[i].Start, r[i].End
	startj, endj := r[j].Start, r[j].End
	if starti != startj {
		return starti < startj
	}
	return endi < endj
}

func (r Refgenes) Swap(i, j int) {
	r[i], r[j] = r[j], r[i]
}

func ReadRefgeneFile(refgeneFile string) Refgenes {
	refgenes := make(Refgenes, 0)
	if fp, err := os.Open(refgeneFile); err == nil {
		defer func(fp *os.File) {
			err := fp.Close()
			if err != nil {
				log.Panic(err)
			}
		}(fp)
		reader := bufio.NewReader(fp)
		for {
			if line, err := reader.ReadBytes('\n'); err == nil {
				refgene := NewRefgene(string(line))
				if refgene.Chrom == "M" || len(refgene.Chrom) > 2 {
					continue
				}
				refgenes = append(refgenes, refgene)
			} else {
				if err == io.EOF {
					break
				} else {
					log.Panic(err)
				}
			}
		}
	} else {
		log.Panic(err)
	}
	sort.Sort(refgenes)
	return refgenes
}
