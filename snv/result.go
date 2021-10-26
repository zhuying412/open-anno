package snv

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
	for i, j := 0, 0; i < len(inputs) && j < len(refIndexes); {
		snvPos1, snvPos2 := inputs[i].Snv.Range()
		refPos1, refPos2 := refIndexes[j].Range()
		if snvPos1 > refPos2 {
			j++
		} else {
			result := Result{Input: inputs[i]}
			if snvPos2 < refPos1 {
				result.Annotations = NewAnnotationsInIntergeic()
			} else {
				var refgenes gene.Refgenes
				refgenes = refgeneMap.FindMany(refIndexes[j].Transcripts)
				result.Annotations = NewAnnotationsInGene(inputs[i].Snv, refgenes)
				if len(result.Annotations) == 0 {
					result.Annotations = NewAnnotationsInUpDownStream(inputs[i].Snv, refgenes)
				}
				if len(result.Annotations) == 0 {
					result.Annotations = NewAnnotationsInIntergeic()
				}
			}
			results = append(results, result)
			i++
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
