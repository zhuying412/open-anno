package io

import (
	"open-anno/pkg"
	"open-anno/pkg/scheme"
	"sort"
	"strings"

	"github.com/brentp/vcfgo"
)

func NewVCFInfo(sample *vcfgo.SampleGenotype, altIndex int) (scheme.VCFInfo, error) {
	altDepths, err := sample.AltDepths()
	if err != nil {
		return scheme.VCFInfo{}, err
	}
	return scheme.VCFInfo{
		Depth: sample.DP,
		VAF:   float64(altDepths[altIndex]) / float64(sample.DP),
		GQ:    float64(sample.GQ),
	}, nil
}

func NewBaseVCFVariant(row *vcfgo.Variant, altIndex int) scheme.VCFVariant {
	return scheme.VCFVariant{
		Chrom:  pkg.FormatChrom(row.Chrom()),
		Pos:    int(row.Pos),
		Ref:    row.Ref(),
		Alt:    row.Alt()[altIndex],
		Qual:   float64(row.Quality),
		Filter: strings.ReplaceAll(row.Filter, ";", ","),
	}
}

func ReadVCF(infile string) (scheme.Variants, error) {
	var variants scheme.Variants
	reader, err := NewIoReader(infile)
	if err != nil {
		return variants, err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return variants, err
	}
	defer vcfReader.Close()
	for {
		row := vcfReader.Read()
		if row == nil {
			break
		}
		sample := row.Samples[0]
		for i := range row.Alt() {
			vcfVariant := NewBaseVCFVariant(row, i)
			vcfVariant.Info, err = NewVCFInfo(sample, i)
			if err != nil {
				return scheme.Variants{}, err
			}
			variants = append(variants, vcfVariant.Variant())
		}
	}
	sort.Sort(variants)
	return variants, err
}

func ReadTriosVCF(infile, proband, mother, father string) (scheme.Variants, error) {
	var variants scheme.Variants
	reader, err := NewIoReader(infile)
	if err != nil {
		return variants, err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return variants, err
	}
	defer vcfReader.Close()
	samples := vcfReader.Header.SampleNames
	pIndex, mIndex, fIndex := pkg.FindArr(samples, proband), pkg.FindArr(samples, mother), pkg.FindArr(samples, father)
	for {
		row := vcfReader.Read()
		if row == nil {
			break
		}
		for i := range row.Alt() {
			vcfVariant := scheme.TriosVCFVariant{VCFVariant: NewBaseVCFVariant(row, i)}
			if pIndex >= 0 {
				proband := row.Samples[pIndex]
				if pkg.Sum(proband.GT...) <= 0 {
					continue
				}
				vcfVariant.Info, err = NewVCFInfo(proband, i)
				if err != nil {
					return scheme.Variants{}, err
				}
			}
			if mIndex >= 0 {
				mother := row.Samples[mIndex]
				if pkg.Sum(mother.GT...) > 0 {
					vcfVariant.MInfo, err = NewVCFInfo(mother, i)
					if err != nil {
						return scheme.Variants{}, err
					}
				}
			}
			if fIndex >= 0 {
				father := row.Samples[fIndex]
				if pkg.Sum(father.GT...) > 0 {
					vcfVariant.FInfo, err = NewVCFInfo(father, i)
					if err != nil {
						return scheme.Variants{}, err
					}
				}
			}
			variants = append(variants, vcfVariant.Variant())
		}
	}
	sort.Sort(variants)
	return variants, err
}

func ReadMcVCF(infile, proband, mother string) (scheme.Variants, error) {
	var variants scheme.Variants
	reader, err := NewIoReader(infile)
	if err != nil {
		return variants, err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return variants, err
	}
	defer vcfReader.Close()
	samples := vcfReader.Header.SampleNames
	pIndex, mIndex := pkg.FindArr(samples, proband), pkg.FindArr(samples, mother)
	for {
		row := vcfReader.Read()
		if row == nil {
			break
		}
		for i := range row.Alt() {
			vcfVariant := scheme.McVCFVariant{VCFVariant: NewBaseVCFVariant(row, i)}
			if pIndex >= 0 {
				proband := row.Samples[pIndex]
				if pkg.Sum(proband.GT...) <= 0 {
					continue
				}
				vcfVariant.Info, err = NewVCFInfo(proband, i)
				if err != nil {
					return scheme.Variants{}, err
				}
			}
			if mIndex >= 0 {
				mother := row.Samples[mIndex]
				if pkg.Sum(mother.GT...) > 0 {
					vcfVariant.MInfo, err = NewVCFInfo(mother, i)
					if err != nil {
						return scheme.Variants{}, err
					}
				}
			}
			variants = append(variants, vcfVariant.Variant())
		}
	}
	sort.Sort(variants)
	return variants, err
}
