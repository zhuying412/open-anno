package input

import (
	"bufio"
	"fmt"
	"grandanno/db"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

type Cnv struct {
	Chrom      string `json:"chrom"`
	Start      int    `json:"start"`
	End        int    `json:"end"`
	Ref        string `json:"ref"`
	Alt        string `json:"alt"`
	CopyNumber int    `json:"copy_number"`
}

func (c Cnv) SN() string {
	return fmt.Sprintf("%s:%d:%d:%s:%s", c.Chrom, c.Start, c.End, c.Ref, c.Alt)
}

func NewCnv(chrom string, start int, end int, copyNumber int) Cnv {
	var alt string
	if copyNumber > 1 {
		alt = "DUP"
	} else if copyNumber < 1 {
		alt = "DEL"
	} else {
		alt = "DIP"
	}
	return Cnv{
		Chrom:      chrom,
		Start:      start,
		End:        end,
		Ref:        "DIP",
		Alt:        alt,
		CopyNumber: copyNumber,
	}
}

type Cnvs []Cnv

func (c Cnvs) Len() int {
	return len(c)
}

func (c Cnvs) Less(i, j int) bool {
	orderChromi, _ := db.ChromArray.GetByName(c[i].Chrom)
	orderChromj, _ := db.ChromArray.GetByName(c[j].Chrom)
	if orderChromi != orderChromj {
		return orderChromi < orderChromj
	}
	starti, endi := c[i].Start, c[i].End
	startj, endj := c[j].Start, c[j].End
	if starti != startj {
		return starti < startj
	}
	return endi < endj
}

func (c Cnvs) Swap(i, j int) {
	c[i], c[j] = c[j], c[i]
}

func (c Cnvs) FilterByChrom(chrom string) Cnvs {
	cnvs := make(Cnvs, 0)
	for _, cnv := range c {
		if cnv.Chrom == chrom {
			cnvs = append(cnvs, cnv)
		}
	}
	return cnvs
}

func ReadCnvInputFile(inputFile string) (cnvs Cnvs, infoMap OtherInfoMap) {
	if fp, err := os.Open(inputFile); err == nil {
		defer func(fp *os.File) {
			err := fp.Close()
			if err != nil {
				log.Panic(err)
			}
		}(fp)
		reader := bufio.NewReader(fp)
		for {
			if line, err := reader.ReadString('\n'); err == nil {
				line = strings.TrimSpace(line)
				if len(line) == 0 || line[0] == '#' {
					continue
				}
				fields := strings.Split(line, "\t")
				start, err := strconv.Atoi(fields[1])
				if err != nil {
					log.Panic(err)
				}
				end, err := strconv.Atoi(fields[2])
				if err != nil {
					log.Panic(err)
				}
				copyNumber, err := strconv.Atoi(fields[3])
				if err != nil {
					log.Panic(err)
				}
				cnv := NewCnv(fields[0], start, end, copyNumber)
				cnvs = append(cnvs, cnv)
				infoMap[cnv.SN()] = NewOtherInfo(fields[4])

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
	return cnvs, infoMap
}
