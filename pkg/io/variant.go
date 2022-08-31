package io

import (
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/scheme"
	"strconv"
	"strings"

	"github.com/brentp/faidx"
)

type VarScanner struct {
	Scanner[scheme.DBVar]
}

func NewVarScanner(reader Reader) VarScanner {
	scanner := NewScanner[scheme.DBVar](reader)
	return VarScanner{Scanner: scanner}
}

func (this VarScanner) Row() (scheme.Variant, error) {
	fields := strings.Split(this.Text(), "\t")
	variant := scheme.Variant{
		Chrom:     pkg.FormatChrom(fields[0]),
		Ref:       fields[3],
		Alt:       fields[4],
		Otherinfo: fields[5],
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

func ReadVariantMap(annoInputFile string) (map[string]scheme.Variants, error) {
	variants := make(map[string]scheme.Variants)
	reader, err := NewIoReader(annoInputFile)
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
			variants[row.Chrom] = scheme.Variants{row}
		}
	}
	return variants, err
}

func ReadVariants(annoInputFile, vcfFile string) (scheme.Variants, error) {
	variants := make(scheme.Variants, 0)
	reader, err := NewIoReader(annoInputFile)
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
		variants = append(variants, row)
	}
	return variants, err
}

func WriteVariants(variants scheme.Variants, outfile string) error {
	writer, err := NewIoWriter(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	for _, row := range variants {
		fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\n", row.Chrom, row.Start, row.End, row.Ref, row.Alt, row.Otherinfo)
	}
	return err
}

func WriteVariantsToVCF(variants scheme.Variants, fai *faidx.Faidx, outfile string) error {
	writer, err := NewIoWriter(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprint(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for _, row := range variants {
		vcfText, err := row.VCFText(fai)
		if err != nil {
			return err
		}
		fmt.Fprintf(writer, vcfText)
	}
	return err
}
