package db

import (
	"bufio"
	"fmt"
	"grandanno/seq"
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
		log.Panic(err.Error())
	}
	if end, err := strconv.Atoi(field[5]); err == nil {
		refgene.End = end
	} else {
		log.Panic(err.Error())
	}
	refgene.UpStream = refgene.Start - UpDownStreamLen
	refgene.DownStream = refgene.End + UpDownStreamLen
	return refgene
}

func (r Refgene) SN() string {
	return fmt.Sprintf("%s|%s:%d:%d", r.Transcript, r.Chrom, r.Start, r.End)
}

func (r Refgene) Range() (int, int) {
	order, _ := ChromArray.GetByName(r.Chrom)
	start := order*1e9 + r.UpStream
	end := order*1e9 + r.DownStream
	return start, end
}

type Refgenes []Refgene

func (refgenes Refgenes) Len() int {
	return len(refgenes)
}

func (refgenes Refgenes) Less(i, j int) bool {
	starti, endi := refgenes[i].Range()
	startj, endj := refgenes[j].Range()
	if starti == startj {
		return endi < endj
	} else {
		return starti < startj
	}
}

func (refgenes Refgenes) Swap(i, j int) {
	refgenes[i], refgenes[j] = refgenes[j], refgenes[i]
}

func NewRefgenes(refgeneFile string) Refgenes {
	refgenes := make(Refgenes, 0)
	if fp, err := os.Open(refgeneFile); err == nil {
		defer func(fp *os.File) {
			err := fp.Close()
			if err != nil {
				log.Panic(err.Error())
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
					log.Panic(err.Error())
				}
			}
		}
	} else {
		log.Panic(err.Error())
	}
	sort.Sort(refgenes)
	return refgenes
}

func WriteMranFasta(refgenes Refgenes, reference seq.Fasta, mrnaFile string) {
	if fp, err := os.Create(mrnaFile); err == nil {
		defer func(fp *os.File) {
			err := fp.Close()
			if err != nil {
				log.Panic(err.Error())
			}
		}(fp)
		for _, refgene := range refgenes {
			if chromSeq, ok := reference[refgene.Chrom]; ok {
				if _, err := fp.WriteString(">" + refgene.SN() + "\n"); err != nil {
					log.Panic(err.Error())
				}
				if _, err := fp.WriteString(string(chromSeq.SubSeq(refgene.Start-1, refgene.End-refgene.Start+1)) + "\n"); err != nil {
					log.Panic(err.Error())
				}
			}
		}
	}
}
