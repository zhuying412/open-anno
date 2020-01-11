package prepare

import (
	"bufio"
	"bytes"
	"fmt"
	"grandanno/core"
	"io"
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

func (refgenes Refgenes) Write(chromSeq []byte, outHandle *os.File) {
	sort.Sort(refgenes)
	for _, refgene := range refgenes {
		refgene.SetSequence(chromSeq)
		if _, err := outHandle.WriteString(">" + refgene.GetSn() + "\n"); err != nil {
			panic(err.Error())
		}
		if _, err := outHandle.WriteString(string(refgene.Sequence) + "\n"); err != nil {
			panic(err.Error())
		}
	}
}

func (refgeneDict RefgeneDict) Read(refgeneFile string, upDownSteamLen int) {
	if fp, err := os.Open(refgeneFile); err == nil {
		defer fp.Close()
		reader := bufio.NewReader(fp)
		for {
			if line, err := reader.ReadBytes('\n'); err == nil {
				var refgene Refgene
				if refgene.Chrom == "M" || len(refgene.Chrom) > 2 {
					continue
				}
				refgene.Read(string(line), upDownSteamLen)
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

func (refgeneDict RefgeneDict) Write(referenceFile string, mrnaFile string) {
	fo, erro := os.Create(mrnaFile)
	fi, erri := os.Open(referenceFile)
	reader := bufio.NewReader(fi)
	if erro == nil && erri == nil {
		defer fo.Close()
		defer fi.Close()
		var name, seq bytes.Buffer
		for {
			if line, err := reader.ReadBytes('\n'); err == nil {
				line = bytes.TrimSpace(line)
				if len(line) == 0 {
					continue
				}
				if line[0] == '>' {
					if name.Len() != 0 {
						chrom := strings.Split(name.String(), " ")[0]
						if refgenes, ok := refgeneDict[chrom]; ok {
							refgenes.Write(seq.Bytes(), fo)
						}
					}
					name.Reset()
					seq.Reset()
					name.Write(line[1:])
				} else {
					seq.Write(line)
				}
			} else {
				if err == io.EOF {
					break
				} else {
					panic(err.Error())
				}
			}
		}
		chrom := strings.Split(name.String(), " ")[0]
		if refgenes, ok := refgeneDict[chrom]; ok {
			refgenes.Write(seq.Bytes(), fo)
		}
	}
}
