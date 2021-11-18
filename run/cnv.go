package run

import (
	"OpenAnno/pkg/variant"
	"bufio"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

func ReadCnvFile(cnvFile string) (variant.Cnvs, OtherInfoMap) {
	cnvs := make(variant.Cnvs, 0)
	infoMap := make(OtherInfoMap)
	fi, err := os.Open(cnvFile)
	if err == nil {
		log.Panic(err)
	}
	defer func(fp *os.File) {
		err := fp.Close()
		if err != nil {
			log.Panic(err)
		}
	}(fi)
	reader := bufio.NewReader(fi)
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Panic(err)
			}
		}
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
		cnv := variant.NewCnv(fields[0], start, end, copyNumber)
		cnvs = append(cnvs, cnv)
		infoMap[cnv.SN()] = NewOtherInfo(fields[5])
	}
	sort.Sort(cnvs)
	return cnvs, infoMap
}
