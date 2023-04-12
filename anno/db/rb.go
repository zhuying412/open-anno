package db

import (
	"fmt"
	"open-anno/anno"
	"open-anno/pkg"
	"path"
	"strings"

	"github.com/brentp/bix"
	"github.com/brentp/irelate/interfaces"
	"github.com/brentp/vcfgo"
)

func annoRegionBased(variants []*pkg.Variant, tbx *bix.Bix, overlap float64, dbname string, annoInfosChan chan anno.AnnoInfos, errChan chan error) {
	annoInfos := make(anno.AnnoInfos)
	for _, variant := range variants {
		query, err := tbx.Query(variant)
		if err != nil {
			annoInfosChan <- anno.AnnoInfos{}
			errChan <- err
			return
		}
		infos := make([]string, 0)
		for v, e := query.Next(); e == nil; v, e = query.Next() {
			info := strings.Split(fmt.Sprintf("%s", v), "\t")[3]
			v := v.(interfaces.IPosition)
			if variant.End() >= v.Start() && variant.Start() <= v.End() {
				vlen := variant.End() - variant.Start()
				olen := pkg.Min(variant.End(), v.End()) - pkg.Max(variant.Start(), v.Start())
				if float64(olen)/float64(vlen) >= overlap {
					infos = append(infos, info)
				}
			}
		}
		annoInfos[variant.PK()] = map[string]any{dbname: strings.Join(infos, "|")}
		query.Close()
	}
	annoInfosChan <- annoInfos
	errChan <- nil
	return
}

// AnnoRegion注释SNV FilterBased
func AnnoRegionBased(inVcfFile, dbBedFile string, overlap float64, goroutines int) (anno.AnnoResult, error) {
	annoInfos := make(anno.AnnoInfos)
	// 打开句柄
	reader, err := pkg.NewIOReader(inVcfFile)
	if err != nil {
		return anno.AnnoResult{}, err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return anno.AnnoResult{}, err
	}
	defer vcfReader.Close()
	dbTbx, err := bix.New(dbBedFile)
	if err != nil {
		return anno.AnnoResult{}, err
	}
	defer dbTbx.Close()
	dbname := strings.Split(path.Base(dbBedFile), ".")[0]
	annoInfosChan := make(chan anno.AnnoInfos, goroutines)
	errChan := make(chan error, goroutines)
	variants := make([]*pkg.Variant, 0)
	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
		if len(variant.Chrom()) <= 5 {
			variants = append(variants, &pkg.Variant{Variant: *variant})
		}
	}
	multiVariants := pkg.SplitArr(variants, goroutines)
	for _, variants := range multiVariants {
		go annoRegionBased(variants, dbTbx, overlap, dbname, annoInfosChan, errChan)
	}
	for i := 0; i < len(multiVariants); i++ {
		err := <-errChan
		if err != nil {
			return anno.AnnoResult{}, err
		}
		for pk, annoInfo := range <-annoInfosChan {
			annoInfos[pk] = annoInfo
		}
	}
	close(annoInfosChan)
	close(errChan)
	vcfHeaderInfos := map[string]*vcfgo.Info{
		dbname: {
			Id:          dbname,
			Description: dbname,
			Number:      ".",
			Type:        "String",
		},
	}
	return anno.AnnoResult{AnnoInfos: annoInfos, VcfHeaderInfo: vcfHeaderInfos}, err
}
