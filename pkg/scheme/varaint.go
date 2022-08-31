package scheme

import (
	"errors"
	"fmt"
	"open-anno/pkg/seq"

	"github.com/brentp/faidx"
)

const (
	VType_SNP = "SNP"
	VType_INS = "INS"
	VType_DEL = "DEL"
	VType_DUP = "DUP"
	VType_SUB = "SUB"
)
const (
	VCMP_GT = ">"
	VCMP_LT = "<"
	VCMP_EQ = "="
)

type Variant struct {
	Chrom     string `json:"chrom"`
	Start     int    `json:"start"`
	End       int    `json:"end"`
	Ref       string `json:"ref"`
	Alt       string `json:"alt"`
	Otherinfo string `json:"otherinfo"`
}

func (this Variant) Type() string {
	if this.Ref == "DIP" {
		if this.Alt == "DEL" {
			return VType_DEL
		}
		return VType_DUP
	}
	if this.Ref == "-" {
		return VType_INS
	} else if this.Alt == "-" {
		return VType_DEL
	} else {
		if len(this.Ref) > 1 || len(this.Alt) > 1 {
			return VType_SUB
		}
		return VType_SNP
	}
}

func (this Variant) ID() string {
	return fmt.Sprintf("%s:%d:%d:%s:%s", this.Chrom, this.Start, this.End, this.Ref, this.Alt)
}

func (this Variant) VCFText(fai *faidx.Faidx) (string, error) {
	chrom, start, ref, alt := this.Chrom, this.Start, this.Ref, this.Alt
	var err error
	if ref != "DIP" {
		if ref == "-" && alt == "-" {
			return "", errors.New("ref == '-' and alt == '-'")
		}
		if ref == "-" {
			ref, err = seq.Fetch(fai, chrom, start-1, start)
			if err == nil {
				alt = ref + alt
			}
		} else if alt == "-" {
			start--
			alt, err = seq.Fetch(fai, chrom, start-1, start)
			if err == nil {
				ref = alt + ref
			}
		}
	}
	return fmt.Sprintf("%s\t%d\t%s\t%s\t%s\t.\t.\t%s", chrom, start, this.ID(), ref, alt, this.Otherinfo), nil
}

func (this Variant) Compare(that Variant) string {
	chrom1, start1, end1, ref1, alt1 := this.Chrom, this.Start, this.End, this.Ref, this.Alt
	chrom2, start2, end2, ref2, alt2 := that.Chrom, that.Start, that.End, that.Ref, that.Alt
	if chrom1 == chrom2 {
		if start1 == start2 {
			if end1 == end2 {
				if ref1 == ref2 {
					if alt1 == alt2 {
						return VCMP_EQ
					}
					if alt1 < alt2 {
						return VCMP_LT
					}
					return VCMP_GT
				}
				if ref1 < ref2 {
					return VCMP_LT
				}
				return VCMP_GT
			}
			if end1 < end2 {
				return VCMP_LT
			}
			return VCMP_GT
		}
		if start1 < start2 {
			return VCMP_LT
		}
		return VCMP_GT

	}
	if chrom1 < chrom2 {
		return VCMP_LT
	}
	return VCMP_GT
}

type Variants []Variant

func (this Variants) Len() int           { return len(this) }
func (this Variants) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this Variants) Less(i, j int) bool { return this[i].Compare(this[j]) == VCMP_LT }
