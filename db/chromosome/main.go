package chromosome

import (
	"bufio"
	"bytes"
	"github.com/spf13/viper"
	"io"
	"log"
	"os"
	"path"
	"strconv"
	"strings"
)

var ChromList Chromosomes

func Init() {
	if len(ChromList) == 0 {
		refDictFile := path.Join(viper.GetString("db.path"), viper.GetString("db.reference_dict"))
		log.Printf("read %s\n", refDictFile)
		ChromList = ReadReferenceDictFile(refDictFile)
	}
}

func ReadReferenceDictFile(referenceDictPath string) Chromosomes {
	chroms := make(Chromosomes, 0)
	fi, err := os.Open(referenceDictPath)
	if err != nil {
		log.Panic(err)
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
				log.Panic(err)
			}
		}
	}
	return chroms
}
