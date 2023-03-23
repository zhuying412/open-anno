package anno

import "github.com/brentp/vcfgo"

type AnnoInfos map[string]map[string]any

type AnnoResult struct {
	AnnoInfos     AnnoInfos
	VcfHeaderInfo map[string]*vcfgo.Info
}

func MergeAnnoResults(annoResults []AnnoResult) AnnoResult {
	annoResult := AnnoResult{AnnoInfos: make(AnnoInfos), VcfHeaderInfo: make(map[string]*vcfgo.Info)}
	for i := 0; i < len(annoResults); i++ {
		annoInfos := annoResults[i].AnnoInfos
		for pk, annoInfo := range annoInfos {
			for key, val := range annoInfo {
				if _, ok := annoResult.AnnoInfos[pk]; !ok {
					annoResult.AnnoInfos[pk] = make(map[string]any)
				}
				annoResult.AnnoInfos[pk][key] = val
			}
		}
		vcfHeaderInfo := annoResults[i].VcfHeaderInfo
		for key, info := range vcfHeaderInfo {
			annoResult.VcfHeaderInfo[key] = info
		}

	}
	return annoResult
}
