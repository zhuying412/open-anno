package pkg

import (
	"fmt"
	"strings"

	"github.com/brentp/vcfgo"
)

const (
	VType_SNP = "SNP"
	VType_INS = "INS"
	VType_DEL = "DEL"
	VType_DUP = "DUP"
	VType_SUB = "SUB"
)

type Position struct {
	chrom string
	start int
	end   int
}

func (this Position) Chrom() string {
	return this.chrom
}
func (this Position) Start() uint64 {
	return uint64(this.start)
}
func (this Position) End() uint64 {
	return uint64(this.end)
}

type AnnoVariant struct {
	Chrom string `json:"chrom"`
	Start int    `json:"start"`
	End   int    `json:"end"`
	Ref   string `json:"ref"`
	Alt   string `json:"alt"`
}

type Variant struct {
	vcfgo.Variant
}

func (this *Variant) PK() string {
	return fmt.Sprintf("%s:%d:%s:%s", this.Chrom(), this.Pos, this.Ref(), this.Alt()[0])
}

func (this *Variant) Type() string {
	ref, alt := this.Ref(), this.Alt()[0]
	if ref == "DIP" {
		return alt
	}
	if len(ref) == 1 {
		if len(alt) == 1 {
			return VType_SNP
		}
		return VType_INS
	} else {
		if len(alt) == 1 {
			return VType_DEL
		}
		return VType_SUB
	}
}

func (this *Variant) AnnoVariant() AnnoVariant {
	chrom, start, ref, alt := this.Chrom(), int(this.Pos), strings.ToUpper(this.Ref()), strings.ToUpper(this.Alt()[0])
	if len(ref) > 1 || len(alt) > 1 && ref != alt {
		if strings.HasPrefix(ref, alt) || strings.HasSuffix(ref, alt) {
			if strings.HasPrefix(ref, alt) {
				start += len(alt)
			}
			ref = strings.Replace(ref, alt, "", 1)
			alt = ""
		} else if strings.HasPrefix(alt, ref) || strings.HasSuffix(alt, ref) {
			if strings.HasPrefix(alt, ref) {
				start += len(ref) - 1
			} else {
				start += len(ref) - len(alt)
			}
			alt = strings.Replace(alt, ref, "", 1)
			ref = ""
		} else {
			refRev := Reverse(ref)
			altRev := Reverse(alt)
			var length int
			length = DifferenceSimple(refRev, altRev) - 1
			ref = ref[0 : len(ref)-length]
			alt = alt[0 : len(alt)-length]
			length = DifferenceSimple(ref, alt) - 1
			ref = ref[length:]
			alt = alt[length:]
			start += length
			if length > 0 && len(ref) == 0 {
				start--
			}
		}
	}
	var end int
	if len(ref) == 0 {
		end = start
		ref = "-"
	} else {
		end = start + len(ref) - 1
	}
	if len(alt) == 0 {
		alt = "-"
	}
	return AnnoVariant{Chrom: chrom, Start: start, End: end, Ref: ref, Alt: alt}
}
