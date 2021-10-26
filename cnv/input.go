package cnv

import (
	"bufio"
	"grandanno/snv"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

type Input struct {
	Cnv       Cnv           `json:"cnv"`
	OtherInfo snv.OtherInfo `json:"other_info"`
}

type Inputs []Input

func (i Inputs) FilterByChrom(chrom string) Inputs {
	inputs := make(Inputs, 0)
	for _, input := range i {
		if input.Cnv.Chrom == chrom {
			inputs = append(inputs, input)
		}
	}
	return inputs
}

func ReadInputFile(inputFile string) (inputs Inputs) {
	if fp, err := os.Open(inputFile); err == nil {
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
				start, err := strconv.Atoi(fields[1])
				if err != nil {
					log.Panic(err)
				}
				end, err := strconv.Atoi(fields[2])
				if err != nil {
					log.Panic(err)
				}
				copyNumber, err := strconv.Atoi(fields[3])
				if err != nil {
					log.Panic(err)
				}
				inputs = append(inputs, Input{
					Cnv:       NewCnv(fields[0], start, end, copyNumber),
					OtherInfo: snv.NewOtherInfo(fields[4]),
				})
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
