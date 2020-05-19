package core

import (
	"fmt"
)

type Variant struct {
	Chrom string
	Start int
	End   int
	Ref   Sequence
	Alt   Sequence
}

func (variant Variant) GetSn() string {
	return fmt.Sprintf("%s:%d:%d:%s:%s", variant.Chrom, variant.Start, variant.End, variant.Ref, variant.Alt)
}

func (variant Variant) GetDigitalPosition() (int, int) {
	start := ChromOrderDict[variant.Chrom]*1e9 + variant.Start
	end := ChromOrderDict[variant.Chrom]*1e9 + variant.End
	return start, end
}

func (variant *Variant) ConvertSnv() {
	if variant.Chrom == "M" {
		variant.Chrom = "MT"
	}
	if !variant.Ref.IsEmpty() || !variant.Alt.IsEmpty() && !variant.Ref.IsEqual(variant.Alt) {
		if variant.Ref.IsStartswith(variant.Alt) || variant.Ref.IsEndswith(variant.Alt) {
			if variant.Ref.IsStartswith(variant.Alt) {
				variant.Start += variant.Alt.GetLen()

			}
			variant.Ref.RemoveOne(variant.Alt)
			variant.Alt.Clear()
		} else if variant.Alt.IsStartswith(variant.Ref) || variant.Alt.IsEndswith(variant.Ref) {
			if variant.Alt.IsStartswith(variant.Ref) {
				variant.Start += len(variant.Ref) - 1
			} else {
				variant.Start += len(variant.Ref) - len(variant.Alt)
			}
			variant.Alt.RemoveOne(variant.Ref)
			variant.Ref.Clear()
		} else {
			var refRev, altRev Sequence
			var subLen int
			refRev = variant.Ref
			altRev = variant.Alt
			refRev.Reverse()
			altRev.Reverse()
			for i, subLen := 0, 0; i < variant.Ref.GetLen() && i < variant.Alt.GetLen(); i++ {
				if refRev.GetChar(i) != altRev.GetChar(i) {
					break
				}
				subLen++
			}
			variant.Ref = variant.Ref.GetSeq(0, variant.Ref.GetLen()-subLen)
			variant.Alt = variant.Alt.GetSeq(0, variant.Alt.GetLen()-subLen)
			for i, subLen := 0, 0; i < variant.Ref.GetLen() && i < variant.Alt.GetLen(); i++ {
				if variant.Ref.GetChar(i) != variant.Alt.GetChar(i) {
					break
				}
				subLen++
			}
			variant.Ref = variant.Ref.GetSeq(subLen, -1)
			variant.Alt = variant.Alt.GetSeq(subLen, -1)
			if subLen > 0 && variant.Ref.IsEmpty() {
				variant.Start += subLen - 1
			} else {
				variant.Start += subLen
			}
		}
	}
	if variant.Ref.IsEmpty() {
		variant.End = variant.Start
		variant.Ref = "-"
	} else {
		variant.End = variant.Start + variant.Ref.GetLen() - 1
	}
	if variant.Alt.IsEmpty() {
		variant.Alt = "-"
	}
}
