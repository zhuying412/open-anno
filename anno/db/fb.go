package db

import (
	"open-anno/anno"
	"open-anno/pkg"

	"github.com/brentp/bix"
	"github.com/brentp/irelate/interfaces"
	"github.com/brentp/vcfgo"
)

// AnnoFilterBased 注释SNV FilterBased
func AnnoFilterBased(snvs []anno.SNV, vcf string) (map[string]map[string]any, []vcfgo.Info, error) {
	annoInfos := make(map[string]map[string]any)
	headerInfos := make([]vcfgo.Info, 0)
	reader, err := pkg.NewIOReader(vcf)
	if err != nil {
		return annoInfos, headerInfos, err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	defer vcfReader.Close()
	for _, info := range vcfReader.Header.Infos {
		headerInfos = append(headerInfos, *info)
	}
	tbx, err := bix.New(vcf)
	if err != nil {
		return annoInfos, headerInfos, err
	}
	defer tbx.Close()
	for _, snv := range snvs {
		iter, err := tbx.Query(&snv)
		if err != nil {
			return annoInfos, headerInfos, err
		}
		for val, err := iter.Next(); err != nil; val, err = iter.Next() {
			variant := val.(interfaces.IVariant)
			if variant.Chrom() == snv.Chrom() && variant.Start() == snv.Start() && variant.End() == snv.End() && variant.Ref() == snv.Ref() && variant.Alt()[0] == snv.Alt()[0] {
				for _, key := range variant.Info().Keys() {
					info, err := variant.Info().Get(key)
					if err != nil {
						return annoInfos, headerInfos, err
					}
					annoInfos[snv.AnnoVariant().PK()] = map[string]any{key: info}
				}

			}
		}
	}
	return annoInfos, headerInfos, err
}
