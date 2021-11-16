package anno

import (
	"OpenAnno/variant"
	"bytes"
	"encoding/json"
	"log"
	"os"
)

type Result struct {
	Variant    variant.IVariant  `json:"variant"`
	Annotation map[string]IAnno  `json:"annotation"`
	OtherInfo  variant.OtherInfo `json:"other_info"`
}

func CreateResultJSON(results []Result, resultFile string) {
	fo, err := os.Create(resultFile)
	if err != nil {
		log.Panic(err)
	}
	defer func(fo *os.File) {
		err = fo.Close()
		if err != nil {
			log.Panic(err)
		}
	}(fo)
	for _, result := range results {
		var buffer bytes.Buffer
		jsonEncoder := json.NewEncoder(&buffer)
		jsonEncoder.SetEscapeHTML(false)
		err = jsonEncoder.Encode(result)
		if err != nil {
			log.Panic(err)
		}
		_, err = fo.Write(buffer.Bytes())
		if err != nil {
			log.Panic(err)
		}
	}
}
