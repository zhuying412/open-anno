package utils

import (
	"encoding/json"
	"log"
	"strconv"
)

func StrToInt(from string) (to int) {
	to, err := strconv.Atoi(from)
	if err != nil {
		log.Panic(err)
	}
	return to
}

func FromJSON(from string, to interface{}) {
	err := json.Unmarshal([]byte(from), &to)
	if err != nil {
		log.Panic(err)
	}
}

func ToJSON(from interface{}) (to string) {
	dat, err := json.Marshal(from)
	if err != nil {
		log.Panic(err)
	}
	return string(dat)
}
