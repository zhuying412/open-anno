package prepare

import (
	"bufio"
	"bytes"
	"fmt"
	"grandanno/bio"
	"grandanno/config"
	"io"
	"os"
	"path"
	"strconv"
	"strings"
)

type Refgene struct {
	Chrom      string
	Start      int
	End        int
	Transcript string
	Sequence   []byte
}

func (refgene *Refgene) Read(refgeneLine string) {
	field := strings.Split(refgeneLine, "\t")
	refgene.Transcript = field[1]
	refgene.Chrom = strings.Replace(field[2], "chr", "", -1)
	if start, err := strconv.Atoi(field[4]); err == nil {
		refgene.Start = start + 1
	}
	if end, err := strconv.Atoi(field[4]); err == nil {
		refgene.End = end
	}
}

func (refgene Refgene) GetSn() string {
	return fmt.Sprintf("%s|%s:%d:%d", refgene.Transcript, refgene.Chrom, refgene.Start, refgene.End)
}

func (refgene Refgene) GetUpStream() int {
	return refgene.Start - config.MyConfig.Param.UpDownStream
}

func (refgene Refgene) GetDownStream() int {
	return refgene.End + config.MyConfig.Param.UpDownStream
}

func (refgene Refgene) GetDigitalPosition() (int, int) {
	start := bio.GetChromOrder(refgene.Chrom)*1e9 + refgene.GetUpStream()
	end := bio.GetChromOrder(refgene.Chrom)*1e9 + refgene.GetDownStream()
	return start, end
}

func (refgene *Refgene) SetSequence(chromSeq []byte) {
	refgene.Sequence = chromSeq[refgene.Start-1 : refgene.End]
}

type Refgenes []Refgene

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

type RefgeneDict map[string]Refgenes

func (refgeneDict RefgeneDict) Read(dbPath string) {
	refgeneFile := path.Join(dbPath, config.MyConfig.Database.Refgene)
	if fp, err := os.Open(refgeneFile); err == nil {
		defer fp.Close()
		reader := bufio.NewReader(fp)
		for {
			if line, err := reader.ReadBytes('\n'); err == nil {
				var refgene Refgene
				refgene.Read(string(line))
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
					panic("read refgene file failed!")
				}
			}
		}
	} else {
		panic("open refgene file failed!")
	}
}

func (refgeneDict RefgeneDict) write(dbPath string) {
	mrnaFile := path.Join(dbPath, config.MyConfig.Database.Mrna)
	referenceFile := path.Join(dbPath, config.MyConfig.Database.Reference)
	fo, erro := os.Create(mrnaFile)
	fi, erri := os.Open(referenceFile)
	reader := bufio.NewReader(fi)
	if erro == nil && erri == nil {
		defer fo.Close()
		defer fi.Close()
		var name, seq bytes.Buffer
		for {
			if line, err := reader.ReadBytes('\n'); err == nil {
				if line[0] == '>' {
					if len(name.Bytes()) != 0 {
						chrom := name.String()
						if refgenes, ok := refgeneDict[chrom]; ok {
							for i := 0; i < len(refgenes); i++ {
								refgenes[i].SetSequence(seq.Bytes())
							}
						}
					}
					name.Write(line[1:])
				} else {
					seq.Write(line)
				}
			} else {
				if err == io.EOF {
					break
				} else {
					panic("read refgene file failed!")
				}
			}
		}
	}
}
