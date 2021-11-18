package variant

import (
	"OpenAnno/pkg/seq"
	"fmt"
	"log"
	"strconv"
)

type VarType string

const (
	VarType_SNV    VarType = "SNV"
	VarType_CNV    VarType = "CNV"
	VarType_COMMON VarType = "COMMON"
)

type VarCmp string

const (
	VarCmp_LT  VarCmp = "LT"
	VarCmp_GT  VarCmp = "GT"
	VarCmp_EQ  VarCmp = "EQ"
	VarCmp_EQP VarCmp = "EQP"
	VarCmp_OV  VarCmp = "OV"
)

type IVariant interface {
	SN() string
	VarType() VarType
	Compare(start int, end int, ref seq.Sequence, alt seq.Sequence) VarCmp
}

type IVariants interface {
	GetVariant(i int) IVariant
	Len() int
	Less(i, j int) bool
	Swap(i, j int)
}

type Variant struct {
	Chrom string       `json:"chrom"`
	Start int          `json:"start"`
	End   int          `json:"end"`
	Ref   seq.Sequence `json:"ref"`
	Alt   seq.Sequence `json:"alt"`
}

func (v Variant) SN() string {
	return fmt.Sprintf("%s:%d:%d:%s:%s", v.Chrom, v.Start, v.End, v.Ref, v.Alt)
}

func (v Variant) VarType() VarType {
	return VarType_COMMON
}

func (v Variant) Compare(start int, end int, ref seq.Sequence, alt seq.Sequence) VarCmp {
	if v.End < start {
		return VarCmp_LT
	} else if v.Start > end {
		return VarCmp_GT
	} else {
		if v.Start == start && v.End == end {
			if v.Ref.IsEqual(ref) && v.Alt.IsEqual(alt) {
				return VarCmp_EQ
			}
			return VarCmp_EQP
		}
		return VarCmp_OV
	}
}

func NewVariant(chrom string, start string, end string, ref string, alt string) Variant {
	_variant := Variant{Chrom: chrom, Ref: seq.Sequence(ref), Alt: seq.Sequence(alt)}
	var err error
	_variant.Start, err = strconv.Atoi(start)
	if err != nil {
		log.Panic(err)
	}
	_variant.End, err = strconv.Atoi(end)
	if err != nil {
		log.Panic(err)
	}
	return _variant
}
