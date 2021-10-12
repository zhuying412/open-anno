package db

import (
	"bufio"
	"bytes"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

type Chromosome struct {
	Name   string
	Length int
}

type Chromosomes []Chromosome

func (c Chromosomes) GetByName(name string) (int, Chromosome) {
	for order, chrom := range c {
		if chrom.Name == name {
			return order, chrom
		}
	}
	log.Panicf("Not found chromosome:%s", name)
	return 0, Chromosome{}
}

func NewChromosomesFromReferenceDict(referenceDictPath string) Chromosomes {
	chroms := make(Chromosomes, 0)
	fi, err := os.Open(referenceDictPath)
	if err != nil {
		log.Panic(err.Error())
	}
	defer func(fi *os.File) {
		err := fi.Close()
		if err != nil {

		}
	}(fi)
	reader := bufio.NewReader(fi)
	for {
		if line, err := reader.ReadBytes('\n'); err == nil {
			line = bytes.TrimSpace(line)
			if len(line) == 0 {
				continue
			}
			fields := bytes.Split(line, []byte{'\t'})
			if bytes.Equal(fields[0], []byte("@SQ")) {
				name := strings.Split(string(fields[1]), ":")[1]
				length, err1 := strconv.Atoi(strings.Split(string(fields[2]), ":")[1])
				if err1 != nil {
					log.Panic(err1.Error())
				}
				chroms = append(chroms, Chromosome{Name: name, Length: length})
			}
		} else {
			if err == io.EOF {
				break
			} else {
				log.Panic(err.Error())
			}
		}
	}
	return chroms
}
