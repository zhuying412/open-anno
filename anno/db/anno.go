package db

import (
	"log"
	"open-anno/anno/variant"
	"open-anno/pkg"
	"strings"
)

// AnnoRegionBased 注释RegionBased
func AnnoSnvs(snvs variant.SNVs, dbFile string, overlap float64) (map[string]map[string]string, error) {
	anno := make(map[string]map[string]string)
	// 读取DBFile DB文件
	log.Printf("Read DBFile: %s ...", dbFile)
	regMap, name, err := ReadRegionBased(dbFile)
	if err != nil {
		return anno, err
	}
	// 开始注释
	for _, snv := range snvs {
		variant := snv.AnnoVariant
		if regs, ok := regMap[variant.Chrom]; ok {
			var infos []string
			for _, reg := range regs {
				if variant.End >= reg.Start && variant.Start <= reg.End {
					vlen := variant.End - variant.Start + 1
					olen := pkg.Min(variant.End, reg.End) - pkg.Max(variant.Start, reg.Start) + 1
					if float64(olen)/float64(vlen) >= overlap {
						infos = append(infos, reg.Info)
					}
				}
			}
			if len(infos) > 0 {
				anno[variant.PK()] = map[string]string{name: strings.Join(infos, ",")}
			}
		}
	}
	return anno, err
}
