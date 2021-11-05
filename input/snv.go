package input

import (
	"bufio"
	"fmt"
	"grandanno/db"
	"grandanno/seq"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

type Snv struct {
	Chrom string       `json:"chrom"`
	Start int          `json:"start"`
	End   int          `json:"end"`
	Ref   seq.Sequence `json:"ref"`
	Alt   seq.Sequence `json:"alt"`
}

func (s Snv) SN() string {
	return fmt.Sprintf("%s:%d:%d:%s:%s", s.Chrom, s.Start, s.End, s.Ref, s.Alt)
}

func (s Snv) Type() string {
	if s.Ref.IsEqual("-") {
		return "ins"
	}
	if s.Alt.IsEqual("-") {
		return "del"
	}
	return "snp"
}

func NewSnv(chrom string, pos int, ref seq.Sequence, alt seq.Sequence) Snv {
	if !ref.IsEmpty() || !alt.IsEmpty() && !ref.IsEqual(alt) {
		if ref.Startswith(alt) || ref.Endswith(alt) {
			if ref.Startswith(alt) {
				pos += alt.Len()
			}
			ref.Replace(alt, 1)
			alt.Clear()
		} else if alt.Startswith(ref) || alt.Endswith(ref) {
			if alt.Startswith(ref) {
				pos += ref.Len() - 1
			} else {
				pos += ref.Len() - alt.Len()
			}
			alt.Replace(ref, 1)
			ref.Clear()
		} else {
			var refRev, altRev seq.Sequence
			var subLen int
			refRev, altRev = ref, alt
			refRev.Reverse()
			altRev.Reverse()
			for i, subLen := 0, 0; i < ref.Len() && i < alt.Len(); i++ {
				if refRev.Base(i) != altRev.Base(i) {
					break
				}
				subLen++
			}
			ref = ref.SubSeq(0, ref.Len()-subLen)
			alt = alt.SubSeq(0, alt.Len()-subLen)
			for i, subLen := 0, 0; i < ref.Len() && i < alt.Len(); i++ {
				if ref.Base(i) != alt.Base(i) {
					break
				}
				subLen++
			}
			ref = ref.SubSeq(subLen, -1)
			alt = alt.SubSeq(subLen, -1)
			if subLen > 0 && ref.IsEmpty() {
				pos += subLen - 1
			} else {
				pos += subLen
			}
		}
	}
	snv := Snv{Chrom: chrom, Start: pos, End: pos, Ref: ref, Alt: alt}
	snv.Chrom = chrom
	if snv.Chrom == "M" {
		snv.Chrom = "MT"
	}
	if snv.Ref.IsEmpty() {
		snv.End = snv.Start
		snv.Ref = "-"
	} else {
		snv.End = snv.Start + snv.Ref.Len() - 1
	}
	if snv.Alt.IsEmpty() {
		snv.Alt = "-"
	}
	return snv
}

type Snvs []Snv

func (s Snvs) Len() int {
	return len(s)
}

func (s Snvs) Less(i, j int) bool {
	orderChromi, _ := db.ChromArray.GetByName(s[i].Chrom)
	orderChromj, _ := db.ChromArray.GetByName(s[j].Chrom)
	if orderChromi != orderChromj {
		return orderChromi < orderChromj
	}
	starti, endi := s[i].Start, s[i].End
	startj, endj := s[j].Start, s[j].End
	if starti != startj {
		return starti < startj
	}
	return endi < endj
}

func (s Snvs) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

type SnvInput struct {
	Snv       Snv       `json:"snv"`
	OtherInfo OtherInfo `json:"other_info"`
}

func (s Snvs) FilterByChrom(chrom string) Snvs {
	snvs := make(Snvs, 0)
	for _, snv := range s {
		if snv.Chrom == chrom {
			snvs = append(snvs, snv)
		}
	}
	return snvs
}

func ReadSnvInputFile(avinputFile string) (Snvs, OtherInfoMap) {
	snvs := make(Snvs, 0)
	infoMap := make(OtherInfoMap)
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
				fields := strings.Split(line, "\t")
				for _, alt := range strings.Split(fields[4], ",") {
					if pos, err := strconv.Atoi(fields[1]); err != nil {
						log.Panic(err)
					} else {
						snv := NewSnv(fields[0], pos, seq.Sequence(fields[3]), seq.Sequence(alt))
						snvs = append(snvs, snv)
						infoMap[snv.SN()] = NewOtherInfo(fields[5])
					}
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
	sort.Sort(snvs)
	return snvs, infoMap
}
