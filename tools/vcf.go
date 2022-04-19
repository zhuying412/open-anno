package tools

import (
	"bufio"
	"errors"
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/seq"
	"open-anno/pkg/variant"
	"os"
	"sort"
	"strings"

	"github.com/brentp/faidx"
	"github.com/brentp/vcfgo"
)

type VCF struct {
	Chrom  string `json:"CHROM"`
	Pos    int    `json:"POS"`
	ID     string `json:"ID"`
	Ref    string `json:"REF"`
	Alt    string `json:"ALT"`
	Qual   string `json:"QUAL"`
	Filter string `json:"FILTER"`
	Info   string `json:"INFO"`
}

type VCFs []VCF

func (this VCFs) Len() int           { return len(this) }
func (this VCFs) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this VCFs) Less(i, j int) bool { return this[i].Pos < this[j].Pos }

func FormatVCFSnv(chrom string, pos int, ref string, alt string) (string, int, int, string, string) {
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

func FormatAVSnv(chrom string, start int, ref string, alt string, fai *faidx.Faidx) (string, int, string, string, error) {
	var pos int
	var err error
	if ref == "-" && alt == "-" {
		return chrom, start, ref, alt, errors.New("ref == '-' and alt == '-'")
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
	return pkg.FormatChrom(chrom), pos, ref, alt, err
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
			chrom, start, end, ref, alt := FormatVCFSnv(pkg.FormatChrom(vcf.Chromosome), int(vcf.Pos), vcf.Ref(), alt)
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

func ReadSnvAV(infile string, fai *faidx.Faidx) (VCFs, error) {
	var vcfs VCFs
	fi, err := os.Open(infile)
	if err != nil {
		return vcfs, err
	}
	defer fi.Close()
	scanner := bufio.NewScanner(fi)
	for scanner.Scan() {
		line := scanner.Text()
		if line[0] == '#' {
			continue
		}
		variant, err := variant.ReadVariantLine(line)
		if err != nil {
			return vcfs, err
		}
		chrom, pos, ref, alt, err := FormatAVSnv(variant.Chrom, variant.Start, variant.Ref, variant.Alt, fai)
		if err != nil {
			return vcfs, err
		}
		vcfs = append(vcfs, VCF{
			Chrom:  chrom,
			Pos:    pos,
			ID:     variant.ID(),
			Ref:    ref,
			Alt:    alt,
			Filter: ".",
			Qual:   ".",
			Info:   strings.Join(variant.Otherinfo, ";"),
		})
	}
	return vcfs, err
}

func WriteVCF(vcfs VCFs, outfile string) error {
	writer, err := os.Create(outfile)
	if err != nil {
		return err
	}
	fmt.Fprint(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for _, vcf := range vcfs {
		fmt.Fprintf(writer, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", vcf.Chrom, vcf.Pos, vcf.ID, vcf.Ref, vcf.Alt, vcf.Qual, vcf.Filter, vcf.Info)
	}
	return nil
}
