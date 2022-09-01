package io

import (
	"open-anno/pkg"
	"open-anno/pkg/scheme"
	"sort"
	"strings"

	"github.com/brentp/vcfgo"
)

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
		alleles, err := sample.AltDepths()
		if err != nil {
			return variants, err
		}
		for i, alt := range row.Alt() {
			vcfVariant := scheme.VCFVariant{
				Chrom:  pkg.FormatChrom(row.Chromosome),
				Pos:    int(row.Pos),
				Ref:    row.Reference,
				Alt:    alt,
				Qual:   float64(row.Quality),
				Filter: strings.ReplaceAll(row.Filter, ";", ","),
				Info: scheme.VCFInfo{
					Depth: sample.DP,
					VAF:   float64(alleles[i]) / float64(sample.DP),
					GQ:    float64(sample.GQ),
				},
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
		var proband, mother, father *vcfgo.SampleGenotype
		var pAlleles, mAlleles, fAlleles []int
		proband = row.Samples[pIndex]
		if pkg.Sum(proband.GT...) <= 0 {
			continue
		}
		pAlleles, err := proband.AltDepths()
		if err != nil {
			return variants, err
		}
		if mIndex >= 0 {
			mother = row.Samples[mIndex]
			mAlleles, err = mother.AltDepths()
			if err != nil {
				return variants, err
			}
		}
		if fIndex >= 0 {
			father = row.Samples[fIndex]
			fAlleles, err = father.AltDepths()
			if err != nil {
				return variants, err
			}
		}
		for i, alt := range row.Alt() {
			vcfVariant := scheme.TriosVCFVariant{
				VCFVariant: scheme.VCFVariant{
					Chrom:  pkg.FormatChrom(row.Chromosome),
					Pos:    int(row.Pos),
					Ref:    row.Reference,
					Alt:    alt,
					Qual:   float64(row.Quality),
					Filter: strings.ReplaceAll(row.Filter, ";", ","),
					Info: scheme.VCFInfo{
						Depth: proband.DP,
						VAF:   float64(pAlleles[1]) / float64(proband.DP),
						GQ:    float64(proband.GQ),
					},
				},
			}
			if mIndex >= 0 && pkg.Sum(mother.GT...) > 0 {
				vcfVariant.MInfo = scheme.VCFInfo{
					Depth: mother.DP,
					VAF:   float64(mAlleles[i]) / float64(mother.DP),
					GQ:    float64(mother.GQ),
				}

			}
			if pIndex >= 0 && pkg.Sum(mother.GT...) > 0 {
				vcfVariant.FInfo = scheme.VCFInfo{
					Depth: father.DP,
					VAF:   float64(fAlleles[i]) / float64(father.DP),
					GQ:    float64(father.GQ),
				}
			}
			variants = append(variants, vcfVariant.Variant())
		}
	}
	sort.Sort(variants)
	return variants, err
}
