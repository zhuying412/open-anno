package tools

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strings"
)

type DBAnno struct {
	Length int
	Header string
	Data   map[string]string
}

func MergeAnno(outfile string, genebased string, otherbaseds ...string) {
	dbAnnos := make([]DBAnno, 0)
	for _, otherbased := range otherbaseds {
		reader, err := os.Open(otherbased)
		if err != nil {
			log.Fatal(err)
		}
		defer reader.Close()
		scanner := bufio.NewScanner(reader)
		scanner.Scan()
		headers := strings.Split(scanner.Text(), "\t")[5:]
		dbAnno := DBAnno{
			Header: strings.Join(headers, "\t"),
			Length: len(headers),
			Data:   map[string]string{},
		}
		for scanner.Scan() {
			fields := strings.Split(scanner.Text(), "\t")
			sn := strings.Join(fields[0:5], ":")
			dbAnno.Data[sn] = strings.Join(fields[5:], "\t")
		}
		dbAnnos = append(dbAnnos, dbAnno)
	}
	reader, err := os.Open(genebased)
	if err != nil {
		log.Fatal(err)
	}
	defer reader.Close()
	writer, err := os.Create(outfile)
	if err != nil {
		log.Fatal(err)
	}
	defer writer.Close()
	scanner := bufio.NewScanner(reader)
	scanner.Scan()
	fmt.Fprint(writer, scanner.Text())
	for _, dbAnno := range dbAnnos {
		fmt.Fprintf(writer, "\t%s", dbAnno.Header)
	}
	fmt.Fprint(writer, "\n")
	for scanner.Scan() {
		line := scanner.Text()
		sn := strings.Join(strings.Split(line, "\t")[0:5], ":")
		fmt.Fprint(writer, line)
		for _, dbAnno := range dbAnnos {
			if anno, ok := dbAnno.Data[sn]; ok {
				fmt.Fprintf(writer, "\t%s", anno)
			} else {
				tmp := make([]string, dbAnno.Length)
				for i := 0; i < dbAnno.Length; i++ {
					tmp[i] = "."
				}
				fmt.Fprintf(writer, "\t%s", strings.Join(tmp, "\t"))
			}
		}
		fmt.Fprint(writer, "\n")
	}
}
