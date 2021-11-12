package database

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
)

func Generate(databaseFile string, outdir string) {
	handlerMap := make(map[string]*os.File)
	if fp, err := os.Open(databaseFile); err == nil {
		defer func(fp *os.File) {
			err := fp.Close()
			if err != nil {
				log.Panic(err)
			}
		}(fp)
		reader := bufio.NewReader(fp)
		var header string
		for {
			if line, err := reader.ReadString('\n'); err == nil {
				line = strings.TrimSpace(line)
				if len(line) == 0 || line[0] == '#' {
					continue
				}
				if header == "" {
					header = line
				}
				fields := strings.Split(line, "\t")
				chrom := fields[0]
				if _, ok := handlerMap[chrom]; !ok {
					handlerMap[chrom], err = os.Create(fmt.Sprintf("%s/chr%s.txt", outdir, chrom))
					if err != nil {
						log.Panic(err)
					}
					_, err = handlerMap[chrom].WriteString(header + "\n")
					if err != nil {
						log.Panic(err)
					}
				}
				_, err = handlerMap[chrom].WriteString(line + "\n")
				if err != nil {
					log.Panic(err)
				}
			} else {
				if err == io.EOF {
					break
				} else {
					log.Panic(err)
				}
			}
		}
	} else {
		log.Panic(err)
	}
	for _, handler := range handlerMap {
		err := handler.Close()
		if err != nil {
			log.Panic(err)
		}
	}
}
