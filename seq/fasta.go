package seq

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
)

type Fasta map[string]Sequence

func (f *Fasta) Set(name string, seq Sequence) {
	name = strings.Split(name, " ")[0]
	(*f)[name] = seq
}

func NewFasta(path string) (fasta Fasta) {
	fasta = make(Fasta)
	fi, err := os.Open(path)
	if err != nil {
		log.Panic(err.Error())
	}
	defer func(fi *os.File) {
		err := fi.Close()
		if err != nil {
			log.Panic(err.Error())
		}
	}(fi)
	reader := bufio.NewReader(fi)
	var name, sequence bytes.Buffer
	for {
		if line, err := reader.ReadBytes('\n'); err == nil {
			line = bytes.TrimSpace(line)
			if len(line) == 0 {
				continue
			}
			if line[0] == '>' {
				if name.Len() != 0 {
					fasta.Set(name.String(), Sequence(sequence.String()))
				}
				name.Reset()
				sequence.Reset()
				name.Write(line[1:])
			} else {
				sequence.Write(line)
			}

		} else {
			if err == io.EOF {
				break
			} else {
				log.Panic(err.Error())
			}
		}
	}
	if name.Len() != 0 {
		fasta.Set(name.String(), Sequence(sequence.String()))
	}
	return fasta
}

func CreateFastaFile(fasta Fasta, path string) {
	fo, err := os.Create(path)
	if err != nil {
		log.Panic(err.Error())
	}
	defer func(fo *os.File) {
		err := fo.Close()
		if err != nil {
			log.Panic(err.Error())
		}
	}(fo)
	for sn, sequence := range fasta {
		if _, err := fo.WriteString(fmt.Sprintf(">%s\n%s\n", sn, sequence)); err != nil {
			log.Panic(err.Error())
		}
	}
}
