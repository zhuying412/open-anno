package schema

import (
	"fmt"
	"open-anno/pkg/seq"
	"strings"
)

type VCFInfo struct {
	VAF   float64 `json:"MVAF"`
	Depth int     `json:"MDEPTH"`
	GQ    float64 `json:"MGQ"`
}

func (this VCFInfo) Text(role string) string {
	return fmt.Sprintf("%sDEPTH=%d;%sVAF=%f;%sGQ=%f;", role, this.Depth, role, this.VAF, role, this.GQ)
}

type IVCFVariant interface {
	ID() string
	InfoText() string
	Variant() Variant
}

type VCFVariant struct {
	Chrom  string  `json:"CHROM"`
	Pos    int     `json:"POS"`
	Ref    string  `json:"REF"`
	Alt    string  `json:"ALT"`
	Qual   float64 `json:"QUAL"`
	Filter string  `json:"FILTER"`
	Info   VCFInfo `json:"INFO"`
}

func (this VCFVariant) ID() string {
	return fmt.Sprintf("%s:%d:%s:%s", this.Chrom, this.Pos, this.Ref, this.Alt)
}

func (this VCFVariant) InfoText() string {
	return fmt.Sprintf("Raw=%s;Qual=%f;Filter=%s;%s",
		this.ID(),
		this.Qual,
		this.Filter,
		this.Info.Text(""),
	)
}

func (this VCFVariant) Variant() Variant {
	chrom, pos, ref, alt := this.Chrom, this.Pos, this.Ref, this.Alt
	if chrom == "M" {
		chrom = "MT"
	}
	start, ref, alt := pos, strings.ToUpper(ref), strings.ToUpper(alt)
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
			refRev := seq.Reverse(ref)
			altRev := seq.Reverse(alt)
			var length int
			length = seq.DifferenceSimple(refRev, altRev) - 1
			ref = ref[0 : len(ref)-length]
			alt = alt[0 : len(alt)-length]
			length = seq.DifferenceSimple(ref, alt) - 1
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
	return Variant{Chrom: chrom, Start: start, End: end, Ref: ref, Alt: alt, Otherinfo: this.InfoText()}
}

type VCFVariants []IVCFVariant

type McVCFVariant struct {
	VCFVariant
	MInfo VCFInfo `json:"MINFO"`
}

func (this McVCFVariant) InfoText() string {
	return fmt.Sprintf("%s;%s",
		this.VCFVariant.InfoText(),
		this.MInfo.Text("M"),
	)
}

func (this McVCFVariant) Variant() Variant {
	variant := this.Variant()
	variant.Otherinfo = this.InfoText()
	return variant
}

type TriosVCFVariant struct {
	VCFVariant
	MInfo VCFInfo `json:"MINFO"`
	FInfo VCFInfo `json:"FINFO"`
}

func (this TriosVCFVariant) InfoText() string {
	return fmt.Sprintf("%s;%s;%s",
		this.VCFVariant.InfoText(),
		this.MInfo.Text("M"),
		this.FInfo.Text("F"),
	)
}

func (this TriosVCFVariant) Variant() Variant {
	variant := this.Variant()
	variant.Otherinfo = this.InfoText()
	return variant
}
