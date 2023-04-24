package pkg

import (
	"fmt"
	"strings"

	"github.com/brentp/irelate/interfaces"
	"github.com/brentp/vcfgo"
)

const (
	VType_SNP = "SNP"
	VType_INS = "INS"
	VType_DEL = "DEL"
	VType_DUP = "DUP"
	VType_SUB = "SUB"
)

type IVariant interface {
	interfaces.IVariant
	PK() string
	AnnoVariant() AnnoVariant
}

type AnnoVariant struct {
	Chrom string `json:"chrom"`
	Start int    `json:"start"`
	End   int    `json:"end"`
	Ref   string `json:"ref"`
	Alt   string `json:"alt"`
}

type SNV struct {
	vcfgo.Variant
}

func (this *SNV) Type() string {
	ref, alt := this.Ref(), this.Alt()[0]
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

func (this *SNV) PK() string {
	return fmt.Sprintf("%s:%d:%s:%s", this.Chrom(), this.Pos, this.Ref(), this.Alt()[0])
}

func (this *SNV) AnnoVariant() AnnoVariant {
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

type CNV struct {
	vcfgo.Variant
}

func (this *CNV) Type() string {
	alt := this.Alt()[0]
	if alt == "<DEL>" || alt == "DEL" || alt == "deletion" {
		return VType_DEL
	}
	return VType_DUP
}

func (this *CNV) PK() string {
	return fmt.Sprintf("%s:%d:%d:%s", this.Chrom(), this.Pos, this.End(), this.Alt()[0])
}

func (this *CNV) AnnoVariant() AnnoVariant {
	chrom, start, end, ref, alt := this.Chrom(), int(this.Pos), int(this.End()), this.Ref(), this.Alt()[0]
	return AnnoVariant{Chrom: chrom, Start: start, End: end, Ref: ref, Alt: alt}
}
