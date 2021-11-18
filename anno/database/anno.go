package database

import (
	"OpenAnno/anno"
	"OpenAnno/pkg/variant"
	"bufio"
	"io"
	"log"
	"strings"
)

type DBAnno map[string]string

func (f DBAnno) AnnoType() anno.AnnoType {
	return anno.AnnoType_DB
}

func ReadLine(reader *bufio.Reader, headers *[]string, vari *variant.Variant, anno *DBAnno) error {
	var fields []string
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				return err
			} else {
				log.Panic(err)
			}
		}
		line = strings.TrimSpace(line)
		if len(line) == 0 || line[0] == '#' {
			if len(*headers) == 0 {
				*headers = strings.Split(line, "\t")
			}
			continue
		}
		fields = strings.Split(line, "\t")
		*vari = variant.NewVariant(fields[0], fields[1], fields[2], fields[3], fields[4])
		for i := 5; i < len(fields); i++ {
			(*anno)[(*headers)[i]] = fields[i]
		}
		return nil
	}
}