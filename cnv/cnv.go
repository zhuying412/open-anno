package cnv

import (
	"bufio"
	"fmt"
	"grandanno/db"
	"io"
	"log"
	"os"
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

func (c Cnv) Range() (int, int) {
	order, _ := db.ChromArray.GetByName(c.Chrom)
	start := order*1e9 + c.Start
	end := order*1e9 + c.End
	return start, end
}

type Cnvs []Cnv

func (c Cnvs) Len() int {
	return len(c)
}

func (c Cnvs) Less(i, j int) bool {
	starti, endi := c[i].Range()
	startj, endj := c[j].Range()
	if starti == startj {
		return endi < endj
	} else {
		return starti < startj
	}
}

func (c Cnvs) Swap(i, j int) {
	c[i], c[j] = c[j], c[i]
}

func NewCnvs(avinputFile string) Cnvs {
	cnvs := make(Cnvs, 0)
	if fp, err := os.Open(avinputFile); err == nil {
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
	return cnvs
}
