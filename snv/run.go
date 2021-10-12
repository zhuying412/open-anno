package snv

import (
	"bytes"
	"encoding/json"
	"grandanno/db"
	"grandanno/gene"
	"log"
	"os"
)

type OpenAnno struct {
	Snv         Snv         `json:"snv"`
	Annotations Annotations `json:"annotation"`
	OtherInfo   interface{} `json:"other_info"`
}

func (o OpenAnno) Json() string {
	var buffer bytes.Buffer
	jsonEncoder := json.NewEncoder(&buffer)
	jsonEncoder.SetEscapeHTML(false)
	if err := jsonEncoder.Encode(o); err != nil {
		log.Panic(err.Error())
	}
	return buffer.String()
}

func RunAnnotate(snvs Snvs, refgeneMap gene.RefgeneMap, refIndexes db.RefIndexes) []OpenAnno {
	annos := make([]OpenAnno, 0)
	for i, j := 0, 0; i < len(snvs) && j < len(refIndexes); {
		snvPos1, snvPos2 := snvs[i].Range()
		refPos1, refPos2 := refIndexes[j].Range()
		if snvPos1 > refPos2 {
			j++
		} else {
			anno := OpenAnno{Snv: snvs[i]}
			if snvPos2 < refPos1 {
				anno.Annotations = NewAnnotationsInIntergeic()
			} else {
				var refgenes gene.Refgenes
				refgenes = refgeneMap.FindMany(refIndexes[j].Transcripts)
				anno.Annotations = NewAnnotationsInGene(snvs[i], refgenes)
				if len(anno.Annotations) == 0 {
					anno.Annotations = NewAnnotationsInUpDownStream(snvs[i], refgenes)
				}
				if len(anno.Annotations) == 0 {
					anno.Annotations = NewAnnotationsInIntergeic()
				}
			}
			annos = append(annos, anno)
			i++
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
