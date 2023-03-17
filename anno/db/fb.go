package db

import (
	"fmt"
	"log"
	"open-anno/anno/gene"
	"open-anno/pkg"

	"github.com/brentp/bix"
	"github.com/brentp/irelate/interfaces"
	"github.com/brentp/vcfgo"
)

func annoFilterBased(variant *pkg.Variant, tbx *bix.Bix, pkChan chan string, annoInfoChan chan map[string]any, errChan chan error) {
	result := make(map[string]any)
	query, err := tbx.Query(variant)
	if err != nil {
		pkChan <- ""
		annoInfoChan <- map[string]any{}
		errChan <- err
		return
	}
	for v, e := query.Next(); e == nil; v, e = query.Next() {
		v := v.(interfaces.IVariant)
		if variant.Chrom() == v.Chrom() && variant.Start() == v.Start() && variant.End() == v.End() && variant.Ref() == v.Ref() && variant.Alt()[0] == v.Alt()[0] {
			for _, key := range v.Info().Keys() {
				val, err := v.Info().Get(key)
				if err == nil {
					result[key] = val
				}
			}
		}
	}
	pkChan <- variant.PK()
	annoInfoChan <- result
	errChan <- nil
	fmt.Println(variant)
	return
}

// AnnoFilterBased 注释SNV FilterBased
func AnnoFilterBased(inVcfFile, dbVcfFile string, goroutines int) (gene.AnnoInfos, map[string]*vcfgo.Info, error) {
	log.Printf("Annotate %s ...", dbVcfFile)
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
	dbTbx, err := bix.New(dbVcfFile)
	fmt.Println("aaaaa", err)
	if err != nil {
		return annoInfos, map[string]*vcfgo.Info{}, err
	}
	defer dbTbx.Close()
	pkChan := make(chan string, goroutines)
	annoInfoChan := make(chan map[string]any, goroutines)
	errChan := make(chan error, goroutines)
	size := 0
	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
		size++
		go annoFilterBased(&pkg.Variant{Variant: *variant}, dbTbx, pkChan, annoInfoChan, errChan)

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
	headerInfos := dbTbx.VReader.Header.Infos
	return annoInfos, headerInfos, err
}
