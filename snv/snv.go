package snv

import (
	"fmt"
	"grandanno/db"
	"grandanno/seq"
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

func (s Snv) Range() (int, int) {
	order, _ := db.ChromArray.GetByName(s.Chrom)
	start := order*1e9 + s.Start
	end := order*1e9 + s.End
	return start, end
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
	starti, endi := s[i].Range()
	startj, endj := s[j].Range()
	if starti == startj {
		return endi < endj
	} else {
		return starti < startj
	}
}

func (s Snvs) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}
