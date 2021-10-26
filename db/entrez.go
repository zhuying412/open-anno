package db

import (
	"bufio"
	"io"
	"log"
	"os"
	"strings"
)

type EntrezIdMap map[string]string

func ReadSymboToIdFile(symboToIdPath string) EntrezIdMap {
	entrezIdMap := make(EntrezIdMap)
	fi, err := os.Open(symboToIdPath)
	if err != nil {
		log.Panic(err)
	}
	defer func(fi *os.File) {
		err := fi.Close()
		if err != nil {
			log.Panic(err)
		}
	}(fi)
	reader := bufio.NewReader(fi)
	for {
		if line, err := reader.ReadString('\n'); err == nil {
			fields := strings.Split(line, "\t")
			symbol, entrezId := fields[0], fields[1]
			entrezIdMap[symbol] = entrezId
		} else {
			if err == io.EOF {
				break
			} else {
				log.Panic(err)
			}
		}
	}
	return entrezIdMap
}

func (e EntrezIdMap) GetEntrezId(symbol string) string {
	if entrezId, ok := e[symbol]; ok {
		return entrezId
	}
	return ""
}
