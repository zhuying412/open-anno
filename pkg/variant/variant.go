package variant

import (
	"bufio"
	"fmt"
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
	return fmt.Sprintf("%s:%d-%d:%s/%s", this.Chrom, this.Start, this.End, this.Ref, this.Alt)
}

func (this Variant) CMP(v Variant) string {
	if this.Start == v.Start {
		if this.End == v.End {
			if this.Ref == v.Ref {
				if this.Alt == v.Alt {
					return VCMP_EQ
				} else {
					if this.Alt < v.Alt {
						return VCMP_LT
					} else {
						return VCMP_GT
					}
				}
			} else {
				if this.Ref < v.Ref {
					return VCMP_LT
				} else {
					return VCMP_GT
				}
			}
		} else {
			if this.End < v.End {
				return VCMP_LT
			} else {
				return VCMP_GT
			}
		}
	} else {
		if this.Start < v.Start {
			return VCMP_LT
		} else {
			return VCMP_GT
		}
	}
}

type Variants []Variant

func (this Variants) Len() int           { return len(this) }
func (this Variants) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this Variants) Less(i, j int) bool { return this[i].CMP(this[j]) == VCMP_LT }

func ReadVariantLine(line string) (Variant, error) {
	fields := strings.Split(line, "\t")
	variant := Variant{
		Chrom:     fields[0],
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
		if vars, ok := variants[variant.Chrom]; ok {
			variants[variant.Chrom] = append(vars, variant)
		} else {
			variants[variant.Chrom] = Variants{variant}
		}
	}
	return variants, err
}
