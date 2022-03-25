package tools

import (
	"fmt"
	"open-anno/pkg/seq"
	"open-anno/pkg/variant"
	"os"
	"sort"
	"strings"

	"github.com/brentp/vcfgo"
)

func FormatSnv(chrom string, pos int, ref string, alt string) (string, int, int, string, string) {
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
			altRef := seq.Reverse(alt)
			var length int
			length = seq.DifferenceSimple(refRev, altRef) - 1
			ref = ref[0 : len(ref)-length]
			alt = ref[0 : len(alt)-length]
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
	return chrom, start, end, ref, alt
}

func ReadVCF(infile string) (variant.Variants, error) {
	var variants variant.Variants
	fi, err := os.Open(infile)
	if err != nil {
		return variants, err
	}
	defer fi.Close()
	reader, err := vcfgo.NewReader(fi, false)
	if err != nil {
		return variants, err
	}
	defer reader.Close()
	for {
		vcf := reader.Read()
		if vcf == nil {
			break
		}
		sample := vcf.Samples[0]
		depth := sample.DP
		alleles, err := sample.AltDepths()
		if err != nil {
			return variants, err
		}
		for i, alt := range vcf.Alt() {
			chrom, start, end, ref, alt := FormatSnv(strings.Replace(vcf.Chromosome, "chr", "", 1), int(vcf.Pos), vcf.Ref(), alt)
			info := fmt.Sprintf("DEPTH=%d;VAF=%f;", depth, float64(alleles[i+1])/float64(depth))
			variants = append(variants, variant.Variant{
				Chrom:     chrom,
				Start:     start,
				End:       end,
				Ref:       ref,
				Alt:       alt,
				Otherinfo: []string{info},
			})
		}
	}
	sort.Sort(variants)
	return variants, err
}
