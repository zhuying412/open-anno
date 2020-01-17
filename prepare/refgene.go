package prepare

import (
	"bufio"
	"fmt"
	"grandanno/core"
	"io"
	"os"
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
	Sequence   []byte
}

type Refgenes []Refgene

type RefgeneDict map[string]Refgenes

func (refgene *Refgene) Read(refgeneLine string, upDownStreamLen int) {
	field := strings.Split(refgeneLine, "\t")
	refgene.Transcript = field[1]
	refgene.Chrom = strings.Replace(strings.Split(field[2], " ")[0], "chr", "", -1)
	if start, err := strconv.Atoi(field[4]); err == nil {
		refgene.Start = start + 1
	}
	if end, err := strconv.Atoi(field[5]); err == nil {
		refgene.End = end
	}
	refgene.UpStream = refgene.Start - upDownStreamLen
	refgene.DownStream = refgene.End + upDownStreamLen
}

func (refgene Refgene) GetSn() string {
	return fmt.Sprintf("%s|%s:%d:%d", refgene.Transcript, refgene.Chrom, refgene.Start, refgene.End)
}

func (refgene Refgene) GetDigitalPosition() (int, int) {
	start := core.ChromOrderDict[refgene.Chrom]*1e9 + refgene.UpStream
	end := core.ChromOrderDict[refgene.Chrom]*1e9 + refgene.DownStream
	return start, end
}

func (refgene *Refgene) SetSequence(chromSeq []byte) {
	refgene.Sequence = chromSeq[refgene.Start-1 : refgene.End]
}

func (refgenes Refgenes) Len() int {
	return len(refgenes)
}

func (refgenes Refgenes) Less(i, j int) bool {
	starti, endi := refgenes[i].GetDigitalPosition()
	startj, endj := refgenes[j].GetDigitalPosition()
	if starti == startj {
		return endi < endj
	} else {
		return starti < startj
	}
}

func (refgenes Refgenes) Swap(i, j int) {
	refgenes[i], refgenes[j] = refgenes[j], refgenes[i]
}

func (refgeneDict RefgeneDict) Read(refgeneFile string, upDownSteamLen int) {
	if fp, err := os.Open(refgeneFile); err == nil {
		defer fp.Close()
		reader := bufio.NewReader(fp)
		for {
			if line, err := reader.ReadBytes('\n'); err == nil {
				var refgene Refgene
				refgene.Read(string(line), upDownSteamLen)
				if refgene.Chrom == "M" || len(refgene.Chrom) > 2 {
					continue
				}
				if refgenes, ok := refgeneDict[refgene.Chrom]; ok {
					refgeneDict[refgene.Chrom] = append(refgenes, refgene)
				} else {
					refgeneDict[refgene.Chrom] = Refgenes{refgene}
				}
			} else {
				if err == io.EOF {
					break
				} else {
					panic(err.Error())
				}
			}
		}
	} else {
		panic(err.Error())
	}
}

func (refgeneDict RefgeneDict) Write(reference core.Fasta, mrnaFile string) {
	if fp, err := os.Create(mrnaFile); err == nil {
		defer fp.Close()
		for _, chrom := range core.ChromList {
			sequence, ok1 := reference[chrom]
			refgenes, ok2 := refgeneDict[chrom]
			if ok1 && ok2 {
				for _, refgene := range refgenes {
					refgene.SetSequence(sequence)
					if _, err := fp.WriteString(">" + refgene.GetSn() + "\n"); err != nil {
						panic(err.Error())
					}
					if _, err := fp.WriteString(string(refgene.Sequence) + "\n"); err != nil {
						panic(err.Error())
					}
				}
			}
		}
	}
}
