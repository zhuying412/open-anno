package cnv

import (
	"encoding/json"
	"grandanno/db"
	"grandanno/gene"
	"log"
	"os"
)

type OpenAnno struct {
	Cnv         Cnv
	Annotations Annotations
	OtherInfo   interface{}
}

func (o OpenAnno) Json() string {
	data, err := json.Marshal(o)
	if err != nil {
		log.Panic(err.Error())
	}
	return string(data)
}

func RunAnnotate(cnvs Cnvs, refgeneMap gene.RefgeneMap, refIndexes db.RefIndexes) []OpenAnno {
	annos := make([]OpenAnno, 0)
	for i := 0; i < len(cnvs); i++ {
		cnvPos1, cnvPos2 := cnvs[i].Range()
		for j := 0; j < len(refIndexes); j++ {
			refPos1, refPos2 := refIndexes[j].Range()
			if cnvPos2 < refPos1 {
				break
			} else if cnvPos1 > refPos2 {
				continue
			} else {
				anno := OpenAnno{Cnv: cnvs[i]}
				if cnvPos2 < refPos1 {
					anno.Annotations = NewAnnotationsInIntergeic()
				} else {
					refgenes := refgeneMap.FindMany(refIndexes[j].Transcripts)
					anno.Annotations = NewAnnotationsInGene(cnvs[i], refgenes)
					if len(anno.Annotations) == 0 {
						anno.Annotations = NewAnnotationsInUpDownStream(cnvs[i], refgenes)
					}
					if len(anno.Annotations) == 0 {
						anno.Annotations = NewAnnotationsInIntergeic()
					}
				}
				annos = append(annos, anno)
			}
		}
	}
	return annos
}

func WriteAnnotations(annos []OpenAnno, outJsonFile string) {
	if fp, err := os.Create(outJsonFile); err == nil {
		defer func(fp *os.File) {
			err := fp.Close()
			if err != nil {
				log.Panic(err.Error())
			}
		}(fp)
		for _, result := range annos {
			if _, err := fp.WriteString(result.Json() + "\n"); err != nil {
				log.Panic(err.Error())
			}
		}
	}
}
