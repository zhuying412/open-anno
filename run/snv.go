package run

import (
	"OpenAnno/pkg/seq"
	"OpenAnno/pkg/variant"
	"bufio"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

func ReadSnvFile(snvFile string) (variant.Snvs, OtherInfoMap) {
	snvs := make(variant.Snvs, 0)
	infoMap := make(OtherInfoMap)
	fi, err := os.Open(snvFile)
	if err != nil {
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
		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			log.Panic(err)
		}
		snv := variant.NewSnv(fields[0], pos, seq.Sequence(fields[3]), seq.Sequence(fields[4]))
		snvs = append(snvs, snv)
		infoMap[snv.SN()] = NewOtherInfo(fields[5])
	}
	sort.Sort(snvs)
	return snvs, infoMap
}
