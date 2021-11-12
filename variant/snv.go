package variant

import (
	"OpenAnno/db/chromosome"
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
	chromOrderi, _ := chromosome.ChromList.GetByName(s[i].Chrom)
	chromOrderj, _ := chromosome.ChromList.GetByName(s[j].Chrom)
	return chromOrderi < chromOrderj || s[i].Start < s[j].Start || s[j].End < s[j].End

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
