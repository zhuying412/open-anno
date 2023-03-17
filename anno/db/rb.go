package db

import (
	"fmt"
	"open-anno/anno/gene"
	"open-anno/pkg"
	"path"
	"strings"

	"github.com/brentp/bix"
	"github.com/brentp/irelate/interfaces"
	"github.com/brentp/vcfgo"
)

func annoRegionBased(variant *pkg.Variant, tbx *bix.Bix, overlap float64, dbname string, pkChan chan string, annoInfoChan chan map[string]any, errChan chan error) {
	query, err := tbx.Query(variant)
	if err != nil {
		pkChan <- ""
		annoInfoChan <- map[string]any{}
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
	pkChan <- variant.PK()
	annoInfoChan <- map[string]any{dbname: strings.Join(infos, ",")}
	errChan <- nil
	return
}

// AnnoRegion注释SNV FilterBased
func AnnoRegionBased(inVcfFile, dbBedFile string, overlap float64, goroutines int) (gene.AnnoInfos, map[string]*vcfgo.Info, error) {
	annoInfos := make(gene.AnnoInfos)
	// 打开句柄
	reader, err := pkg.NewIOReader(inVcfFile)
	if err != nil {
		return annoInfos, map[string]*vcfgo.Info{}, err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return annoInfos, map[string]*vcfgo.Info{}, err
	}
	defer vcfReader.Close()
	dbTbx, err := bix.New(dbBedFile)
	if err != nil {
		return annoInfos, map[string]*vcfgo.Info{}, err
	}
	defer dbTbx.Close()
	pkChan := make(chan string, goroutines)
	annoInfoChan := make(chan map[string]any, goroutines)
	errChan := make(chan error, goroutines)
	size := 0
	dbname := strings.Split(path.Base(dbBedFile), ".")[0]
	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
		size++
		go annoRegionBased(&pkg.Variant{Variant: *variant}, dbTbx, overlap, dbname, pkChan, annoInfoChan, errChan)

	}
	for i := 0; i < size; i++ {
		err = <-errChan
		if err != nil {
			return annoInfos, map[string]*vcfgo.Info{}, err
		}
		annoInfos[<-pkChan] = <-annoInfoChan
	}
	close(pkChan)
	close(annoInfoChan)
	close(errChan)
	headerInfos := map[string]*vcfgo.Info{
		dbname: {
			Id:          dbname,
			Description: dbname,
			Number:      ".",
			Type:        "String",
		},
	}
	return annoInfos, headerInfos, err
}
