package io

import (
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/schema"
	"strconv"
	"strings"

	"github.com/brentp/faidx"
)

type VarScanner struct {
	GenericsTSVScanner[schema.DBVar]
}

func NewVarScanner(reader Reader) VarScanner {
	scanner := NewGenericsTSVScanner[schema.DBVar](reader)
	return VarScanner{GenericsTSVScanner: scanner}
}

func (this VarScanner) Row() (schema.Variant, error) {
	fields := strings.Split(this.Text(), "\t")
	variant := schema.Variant{
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

func ReadVariantMap(annoInputFile string) (map[string]schema.Variants, error) {
	variants := make(map[string]schema.Variants)
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
			variants[row.Chrom] = schema.Variants{row}
		}
	}
	return variants, err
}

func ReadVariants(annoInputFile string) (schema.Variants, error) {
	variants := make(schema.Variants, 0)
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

func WriteVariants(variants schema.Variants, outfile string) error {
	writer, err := NewIoWriter(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprint(writer, "Chr\tStart\tEnd\tRef\tAlt\tOtherinfo\n")
	for _, row := range variants {
		fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\n", row.Chrom, row.Start, row.End, row.Ref, row.Alt, row.Otherinfo)
	}
	return err
}

func WriteVariantsToVCF(variants schema.Variants, fai *faidx.Faidx, outfile string) error {
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
