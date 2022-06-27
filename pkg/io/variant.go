package io

import (
	"errors"
	"fmt"
	"io"
	"open-anno/pkg"
	"open-anno/pkg/seq"
	"strconv"
	"strings"

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

func (this Variant) Compare(that Variant) string {
	chrom1, start1, end1, ref1, alt1 := this.Chrom, this.Start, this.End, this.Ref, this.Alt
	chrom2, start2, end2, ref2, alt2 := that.Chrom, that.Start, that.End, that.Ref, that.Alt
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

func (this Variant) VCF(fai *faidx.Faidx) (VCF, error) {
	chrom, start, end, ref, alt := this.Chrom, this.Start, this.End, this.Ref, this.Alt
	var pos int
	var err error
	if ref != "DIP" {
		if ref == "-" && alt == "-" {
			return VCF{}, errors.New("ref == '-' and alt == '-'")
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
	return VCF{
		Chrom: chrom, Pos: pos, Ref: ref, Alt: alt,
		ID:   fmt.Sprintf("%s:%d:%d:%s:%s", chrom, start, end, ref, alt),
		Qual: 0, Filter: ".", Info: strings.Join(this.Otherinfo, ","),
	}, nil
}

type VarScanner struct {
	Scanner[Variant]
}

func NewVarScanner(reader io.Reader) VarScanner {
	scanner := NewScanner[Variant](reader)
	return VarScanner{Scanner: scanner}
}

func (this VarScanner) Row() (Variant, error) {
	fields := strings.Split(this.Text(), "\t")
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
	return variant, nil
}

type Variants []Variant

func (this Variants) Len() int           { return len(this) }
func (this Variants) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this Variants) Less(i, j int) bool { return this[i].Compare(this[j]) == VCMP_LT }

func ReadVariantMap(avinput string) (map[string]Variants, error) {
	variants := make(map[string]Variants)
	reader, err := NewIoReader(avinput)
	if err != nil {
		return variants, err
	}
	defer reader.Close()
	scanner := NewVarScanner(reader)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return variants, err
		}
		if rows, ok := variants[row.Chrom]; ok {
			variants[row.Chrom] = append(rows, row)
		} else {
			variants[row.Chrom] = Variants{row}
		}
	}
	return variants, err
}

func WriteVariants(outfile string, variants ...Variants) error {
	writer, err := NewIoWriter(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	for _, rows := range variants {
		for _, row := range rows {
			fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\n", row.Chrom, row.Start, row.End, row.Ref, row.Alt, row.Otherinfo)
		}
	}
	return err
}
