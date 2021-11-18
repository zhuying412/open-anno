package variant

import (
	"OpenAnno/pkg/seq"
)

type SnvType string

const (
	SnvType_SNP SnvType = "SNP"
	SnvType_DEL SnvType = "DEL"
	SnvType_INS SnvType = "INS"
)

type Snv struct {
	Variant
}

func (s Snv) VarType() VarType {
	return VarType_SNV
}

func (s Snv) Type() SnvType {
	if s.Ref.IsEqual("-") {
		return SnvType_INS
	}
	if s.Alt.IsEqual("-") {
		return SnvType_DEL
	}
	return SnvType_SNP
}

type Snvs []Snv

func (s Snvs) GetVariant(i int) IVariant {
	return s[i]
}

func (s Snvs) Len() int {
	return len(s)
}

func (s Snvs) Less(i, j int) bool {
	return s[i].Start < s[j].Start || s[j].End < s[j].End

}

func (s Snvs) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
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
	snv := Snv{Variant: Variant{Chrom: chrom, Start: pos, End: pos, Ref: ref, Alt: alt}}
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
