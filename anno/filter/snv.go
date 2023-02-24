package db

import (
	"log"
	"open-anno/anno/variant"
	"sort"
)

// AnnoSnvs 注释SNV FilterBased
func AnnoSnvs(snvs variant.SNVs, dbFile, dbIdxFile, outfile string) (map[string]map[string]string, error) {
	anno := make(map[string]map[string]string)
	// 读取DB Idx File输入文件
	log.Printf("Read DBIdxFile: %s ...", dbIdxFile)
	idxMap, binSize, err := ReadFilterVarIdx(dbIdxFile)
	if err != nil {
		return anno, err
	}
	// 转换SNV列表
	log.Printf("Aggregate SNVs by Bin...")
	snvMap := snvs.AggregateByBin(binSize)
	// 开始注释
	reader, err := NewFilterVarReader(dbFile)
	if err != nil {
		return anno, err
	}
	defer reader.Close()
	for curbin, subsnvs := range snvMap {
		sort.Sort(subsnvs)
		if idx, ok := idxMap[curbin]; ok {
			dbvars, err := reader.ReadRows(idx)
			if err != nil {
				return anno, err
			}
			for i, j := 0, 0; i < len(subsnvs) && j < len(dbvars); {
				if dbvars[j].Greater(subsnvs[i].AnnoVariant) {
					i++
				} else if dbvars[j].Less(subsnvs[i].AnnoVariant) {
					j++
				} else {
					anno[subsnvs[i].AnnoVariant.PK()] = dbvars[j].Info
					i++
				}
			}
		}
	}
	return anno, err
}
