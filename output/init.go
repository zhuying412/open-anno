package output

import (
	"bytes"
	"encoding/json"
	"log"
	"os"
)

type Outputs interface {
	JsonLines() []string
}

func ConvertToJSON(instance interface{}) string {
	var buffer bytes.Buffer
	jsonEncoder := json.NewEncoder(&buffer)
	jsonEncoder.SetEscapeHTML(false)
	if err := jsonEncoder.Encode(instance); err != nil {
		log.Panic(err)
	}
	return buffer.String()
}

func CreateOutputFile(outputs Outputs, outJsonFile string) {
	if fp, err := os.Create(outJsonFile); err == nil {
		defer func(fp *os.File) {
			err := fp.Close()
			if err != nil {
				log.Panic(err)
			}
		}(fp)
		for _, line := range outputs.JsonLines() {
			if _, err = fp.WriteString(line + "\n"); err != nil {
				log.Panic(err)
			}
		}
	}
}
