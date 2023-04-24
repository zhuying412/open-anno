package db

import (
	"open-anno/pkg"

	"github.com/brentp/bix"
	"github.com/brentp/irelate/interfaces"
)

func AnnoFilterBased(variant pkg.IVariant, tbx *bix.Bix) (map[string]any, error) {
	query, err := tbx.Query(variant)
	if err != nil {
		return map[string]any{}, err
	}
	annoInfo := make(map[string]any)
	for v, e := query.Next(); e == nil; v, e = query.Next() {
		v := v.(interfaces.IVariant)
		if variant.Chrom() == v.Chrom() && variant.Start() == v.Start() && variant.End() == v.End() && variant.Ref() == v.Ref() && variant.Alt()[0] == v.Alt()[0] {
			for _, key := range v.Info().Keys() {
				val, err := v.Info().Get(key)
				if err == nil {
					annoInfo[key] = val
				}
			}
		}
	}
	query.Close()
	return annoInfo, err
}

// // AnnoFilterBased 注释SNV FilterBased
// func AnnoFilterBased(inVcfFile, dbVcfFile string, goroutines int) (anno.AnnoResult, error) {
// 	log.Printf("Annotate %s ...", dbVcfFile)
// 	annoInfos := make(anno.AnnoInfos)
// 	// 打开句柄
// 	reader, err := pkg.NewIOReader(inVcfFile)
// 	if err != nil {
// 		return anno.AnnoResult{}, err
// 	}
// 	defer reader.Close()
// 	vcfReader, err := vcfgo.NewReader(reader, false)
// 	if err != nil {
// 		return anno.AnnoResult{}, err
// 	}
// 	defer vcfReader.Close()
// 	dbTbx, err := bix.New(dbVcfFile)
// 	if err != nil {
// 		return anno.AnnoResult{}, err
// 	}
// 	defer dbTbx.Close()
// 	annoInfosChan := make(chan anno.AnnoInfos, goroutines)
// 	errChan := make(chan error, goroutines)
// 	variants := make([]*pkg.Variant, 0)
// 	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
// 		if len(variant.Chrom()) <= 5 {
// 			variants = append(variants, &pkg.Variant{Variant: *variant})
// 		}
// 	}
// 	multiVariants := pkg.SplitArr(variants, goroutines)
// 	for _, variants := range multiVariants {
// 		go annoFilterBased(variants, dbTbx, annoInfosChan, errChan)
// 	}
// 	for i := 0; i < len(multiVariants); i++ {
// 		err := <-errChan
// 		if err != nil {
// 			return anno.AnnoResult{}, err
// 		}
// 		for pk, annoInfo := range <-annoInfosChan {
// 			annoInfos[pk] = annoInfo
// 		}
// 	}
// 	close(annoInfosChan)
// 	close(errChan)
// 	vcfHeaderInfos := dbTbx.VReader.Header.Infos
// 	return anno.AnnoResult{AnnoInfos: annoInfos, VcfHeaderInfo: vcfHeaderInfos}, err
// }
