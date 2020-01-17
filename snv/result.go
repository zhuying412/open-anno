package snv

import (
	"encoding/json"
	"grandanno/core"
	"os"
)

type Result struct {
	Snv         Snv
	Annotations Annotations
}

type Results []Result

func (result Result) GetJson() string {
	if data, err := json.Marshal(result); err == nil {
		return string(data)
	} else {
		panic(err)
	}
}

func (results *Results) RunAnno(snvs Snvs, refgeneDict core.RefgeneDict, refidxs core.Refidxs, splicingLen int) {
	for i, j := 0, 0; i < len(snvs) && j < len(refidxs); {
		snvPos1, snvPos2 := snvs[i].Variant.GetDigitalPosition()
		refPos1, refPos2 := refidxs[j].GetDigitalPosition()
		if snvPos1 > refPos2 {
			j++
		} else {
			result := Result{Snv: snvs[i]}
			if snvPos2 < refPos1 {
				result.Annotations.AnnoIntergeic()
			} else {
				refgenes := refidxs[j].GetRefgenes(refgeneDict)
				result.Annotations.AnnoGene(snvs[i], refgenes, splicingLen)
				if len(result.Annotations) == 0 {
					result.Annotations.AnnoStream(snvs[i], refgenes)
				}
				if len(result.Annotations) == 0 {
					result.Annotations.AnnoIntergeic()
				}
			}
			*results = append(*results, result)
			i++
		}
	}
}

func (results Results) Write(outJsonFile string) {
	if fp, err := os.Create(outJsonFile); err == nil {
		defer fp.Close()
		for _, result := range results {
			if _, err := fp.WriteString(result.GetJson() + "\n"); err != nil {
				panic(err)
			}
		}
	}
}
