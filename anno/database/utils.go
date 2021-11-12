package database

import (
	"OpenAnno/variant"
	"bufio"
	"io"
	"log"
	"strings"
)

func ReadLine(reader *bufio.Reader, headers *[]string, vari *variant.Variant, anno *IAnno) error {
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
		fields = strings.Split(line, "\t")
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		if len(*headers) == 0 {
			*headers = fields
		} else {
			*vari = variant.NewVariant(fields[0], fields[1], fields[2], fields[3], fields[4])
			for i := 5; i < len(fields); i++ {
				(*anno).Set((*headers)[i], fields[i])
			}
			return nil
		}
	}
}
