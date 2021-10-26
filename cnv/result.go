package cnv

import (
	"bytes"
	"encoding/json"
	"grandanno/db"
	"grandanno/gene"
	"log"
	"os"
)

type Result struct {
	Input
	Annotations Annotations `json:"annotation"`
}

func (r Result) Json() string {
	var buffer bytes.Buffer
	jsonEncoder := json.NewEncoder(&buffer)
	jsonEncoder.SetEscapeHTML(false)
	if err := jsonEncoder.Encode(r); err != nil {
		log.Panic(err)
	}
	return buffer.String()
}

func RunAnnotate(inputs Inputs, refgeneMap gene.RefgeneMap, refIndexes db.RefIndexes) []Result {
	results := make([]Result, 0)
	for i := 0; i < len(inputs); i++ {
		cnvPos1, cnvPos2 := inputs[i].Cnv.Range()
		for j := 0; j < len(refIndexes); j++ {
			refPos1, refPos2 := refIndexes[j].Range()
			if cnvPos2 < refPos1 {
				break
			} else if cnvPos1 > refPos2 {
				continue
			} else {
				result := Result{Input: inputs[i]}
				if cnvPos2 < refPos1 {
					result.Annotations = NewAnnotationsInIntergeic()
				} else {
					refgenes := refgeneMap.FindMany(refIndexes[j].Transcripts)
					result.Annotations = NewAnnotationsInGene(inputs[i].Cnv, refgenes)
					if len(result.Annotations) == 0 {
						result.Annotations = NewAnnotationsInUpDownStream(inputs[i].Cnv, refgenes)
					}
					if len(result.Annotations) == 0 {
						result.Annotations = NewAnnotationsInIntergeic()
					}
				}
				results = append(results, result)
			}
		}
	}
	return results
}

func CreateAnnotationFile(results []Result, outJsonFile string) {
	if fp, err := os.Create(outJsonFile); err == nil {
		defer func(fp *os.File) {
			err := fp.Close()
			if err != nil {
				log.Panic(err)
			}
		}(fp)
		for _, result := range results {
			if _, err := fp.WriteString(result.Json() + "\n"); err != nil {
				log.Panic(err)
			}
		}
	}
}
