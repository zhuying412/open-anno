package variant

import (
	"bufio"
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/gene"
	"os"
	"strconv"
	"strings"
)

const (
	VType_SNP = "SNP"
	VType_INS = "INS"
	VType_DEL = "DEL"
	VType_DUP = "DUP"
)
const (
	VCMP_GT = ">"
	VCMP_LT = "<"
	VCMP_EQ = "="
)

type Variant struct {
	Chrom     string   `json:"chrom"`
	Start     int      `json:"start"`
	End       int      `json:"end"`
	Ref       string   `json:"ref"`
	Alt       string   `json:"alt"`
	Otherinfo []string `json:"otherinfo"`
}

func (this Variant) Type() string {
	if this.Ref == "DIP" {
		if this.Alt == "DEL" {
			return VType_DEL
		}
		return VType_DUP
	}
	if this.Alt == "-" {
		return VType_DEL
	} else if this.Ref == "-" {
		return VType_INS
	}
	return VType_SNP
}

func (this Variant) ID() string {
	return fmt.Sprintf("%s:%d:%d:%s:%s", this.Chrom, this.Start, this.End, this.Ref, this.Alt)
}

func (this Variant) GetBaseVar() (string, int, int, string, string) {
	return this.Chrom, this.Start, this.End, this.Ref, this.Alt
}

type Variants []Variant

func (this Variants) Len() int           { return len(this) }
func (this Variants) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this Variants) Less(i, j int) bool { return CompareVar(this[i], this[j]) == VCMP_LT }

func ReadVariantLine(line string) (Variant, error) {
	fields := strings.Split(line, "\t")
	variant := Variant{
		Chrom:     pkg.FormatChrom(fields[0]),
		Ref:       fields[3],
		Alt:       fields[4],
		Otherinfo: fields[5:],
	}
	var err error
	variant.Start, err = strconv.Atoi(fields[1])
	if err != nil {
		return variant, err
	}
	variant.End, err = strconv.Atoi(fields[2])
	if err != nil {
		return variant, err
	}
	return variant, err
}

func ReadAvinput(avinput string) (map[string]Variants, error) {
	variants := make(map[string]Variants)
	fi, err := os.Open(avinput)
	if err != nil {
		return variants, err
	}
	defer fi.Close()
	scanner := bufio.NewScanner(fi)
	for scanner.Scan() {
		line := scanner.Text()
		if line[0] == '#' {
			continue
		}
		variant, err := ReadVariantLine(line)
		if err != nil {
			return variants, err
		}
		if _, ok := gene.GENOME[variant.Chrom]; !ok {
			continue
		}
		if vars, ok := variants[variant.Chrom]; ok {
			variants[variant.Chrom] = append(vars, variant)
		} else {
			variants[variant.Chrom] = Variants{variant}
		}
	}
	return variants, err
}

func WriteAvinput(variants Variants, outfile string) error {
	writer, err := os.Create(outfile)
	if err != nil {
		return err
	}
	for _, variant := range variants {
		fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s", variant.Chrom, variant.Start, variant.End, variant.Ref, variant.Alt)
		for _, info := range variant.Otherinfo {
			fmt.Fprintf(writer, "\t%s", info)
		}
		fmt.Fprint(writer, "\n")
	}
	return nil
}

type ICompareVar interface {
	GetBaseVar() (string, int, int, string, string)
}

func CompareVar(v1 ICompareVar, v2 ICompareVar) string {
	chrom1, start1, end1, ref1, alt1 := v1.GetBaseVar()
	chrom2, start2, end2, ref2, alt2 := v2.GetBaseVar()
	if chrom1 == chrom2 {
		if start1 == start2 {
			if end1 == end2 {
				if ref1 == ref2 {
					if alt1 == alt2 {
						return VCMP_EQ
					} else {
						if alt1 < alt2 {
							return VCMP_LT
						} else {
							return VCMP_GT
						}
					}
				} else {
					if ref1 < ref2 {
						return VCMP_LT
					} else {
						return VCMP_GT
					}
				}
			} else {
				if end1 < end2 {
					return VCMP_LT
				} else {
					return VCMP_GT
				}
			}
		} else {
			if start1 < start2 {
				return VCMP_LT
			} else {
				return VCMP_GT
			}
		}
	} else {
		if chrom1 < chrom2 {
			return VCMP_LT
		} else {
			return VCMP_GT
		}
	}
}
