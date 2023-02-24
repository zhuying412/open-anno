package variant

import (
	"open-anno/pkg"

	"github.com/brentp/vcfgo"
)

type SNV struct {
	vcfgo.Variant
	AnnoVariant AnnoVariant
}

type SNVs []SNV

func (this SNVs) Len() int           { return len(this) }
func (this SNVs) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this SNVs) Less(i, j int) bool { return this[i].AnnoVariant.Less(this[j].AnnoVariant) }

func (this SNVs) AggregateByChrom() map[string]SNVs {
	snvs := make(map[string]SNVs)
	for _, snv := range this {
		chrom := snv.AnnoVariant.Chrom
		if rows, ok := snvs[chrom]; ok {
			snvs[chrom] = append(rows, snv)
		} else {
			snvs[chrom] = SNVs{snv}
		}
	}
	return snvs
}

func (this SNVs) AggregateByBin(binSize int) map[string]SNVs {
	snvs := make(map[string]SNVs)
	for _, snv := range this {
		chrom, start := snv.AnnoVariant.Chrom, snv.AnnoVariant.Start
		curbin := pkg.CurBin(chrom, start, binSize)
		if rows, ok := snvs[curbin]; ok {
			snvs[curbin] = append(rows, snv)
		} else {
			snvs[curbin] = SNVs{snv}
		}
	}
	return snvs
}

func ReadVCF(vcfFile string) (SNVs, error) {
	snvs := make(SNVs, 0)
	reader, err := pkg.NewIOReader(vcfFile)
	if err != nil {
		return snvs, err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return snvs, err
	}
	defer vcfReader.Close()
	for variant := vcfReader.Read(); variant != nil; {
		for _, alt := range variant.Alt() {
			chrom, start, end, ref, alt := pkg.VCFtoAV(variant.Chrom(), int(variant.Pos), variant.Ref(), alt)
			snv := SNV{
				Variant:     *variant,
				AnnoVariant: AnnoVariant{Chrom: chrom, Start: start, End: end, Ref: ref, Alt: alt},
			}
			snvs = append(snvs, snv)
		}

	}
	return snvs, nil
}
