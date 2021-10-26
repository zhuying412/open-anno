package cnv

import (
	"fmt"
	"grandanno/db"
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
