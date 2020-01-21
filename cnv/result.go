package cnv

import (
	"encoding/json"
	"grandanno/core"
	"os"
)

type Result struct {
	Cnv         Cnv
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

func (results *Results) RunAnno(cnvs Cnvs, refgeneDict core.RefgeneDict, refidxs core.Refidxs) {
	for i := 0; i < len(cnvs); i++ {
		cnvPos1, cnvPos2 := cnvs[i].GetVariant().GetDigitalPosition()
		for j := 0; j < len(refidxs); j++ {
			refPos1, refPos2 := refidxs[j].GetDigitalPosition()
			if cnvPos2 < refPos1 {
				break
			} else if cnvPos1 > refPos2 {
				continue
			} else {
				result := Result{Cnv: cnvs[i]}
				if cnvPos2 < refPos1 {
					result.Annotations.AnnoIntergeic()
				} else {
					refgenes := refidxs[j].GetRefgenes(refgeneDict)
					result.Annotations.AnnoGene(cnvs[i], refgenes)
					if len(result.Annotations) == 0 {
						result.Annotations.AnnoStream(cnvs[i], refgenes)
					}
					if len(result.Annotations) == 0 {
						result.Annotations.AnnoIntergeic()
					}
				}
			}
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
