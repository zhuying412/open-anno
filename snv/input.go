package snv

import (
	"bufio"
	"grandanno/seq"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

type Input struct {
	Snv       Snv       `json:"snv"`
	OtherInfo OtherInfo `json:"other_info"`
}

type Inputs []Input

func (i Inputs) FilterByChrom(chrom string) Inputs {
	inputs := make(Inputs, 0)
	for _, input := range i {
		if input.Snv.Chrom == chrom {
			inputs = append(inputs, input)
		}
	}
	return inputs
}

func ReadInputFile(avinputFile string) (inputs Inputs) {
	if fp, err := os.Open(avinputFile); err == nil {
		defer func(fp *os.File) {
			err := fp.Close()
			if err != nil {
				log.Panic(err)
			}
		}(fp)
		reader := bufio.NewReader(fp)
		for {
			if line, err := reader.ReadString('\n'); err == nil {
				line = strings.TrimSpace(line)
				if len(line) == 0 || line[0] == '#' {
					continue
				}
				fields := strings.Split(line, "\t")
				for _, alt := range strings.Split(fields[4], ",") {
					if pos, err := strconv.Atoi(fields[1]); err != nil {
						log.Panic(err)
					} else {
						inputs = append(inputs,
							Input{
								Snv:       NewSnv(fields[0], pos, seq.Sequence(fields[3]), seq.Sequence(alt)),
								OtherInfo: NewOtherInfo(fields[5]),
							},
						)
					}
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
	return inputs
}
