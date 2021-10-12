package db

import (
	"bufio"
	"fmt"
	"grandanno/seq"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

func NewMrna(refgenePath string, chromSeq seq.Fasta) seq.Fasta {
	mrna := make(seq.Fasta)
	fi, err := os.Open(refgenePath)
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
	for {
		if line, err := reader.ReadString('\n'); err == nil {
			fields := strings.Split(line, "\t")
			transcript := fields[1]
			chrom := strings.Replace(strings.Split(fields[2], " ")[0], "chr", "", -1)
			start, err := strconv.Atoi(fields[4])
			if err != nil {
				log.Panic(err.Error())
			}
			start = start + 1
			end, err := strconv.Atoi(fields[5])
			if err != nil {
				log.Panic(err.Error())
			}
			if chrom == "M" || len(chrom) > 2 {
				continue
			}
			if sequence, ok := chromSeq[chrom]; ok {
				sn := fmt.Sprintf("%s|%s:%d:%d", transcript, chrom, start, end)
				mrna[sn] = sequence.SubSeq(start-1, end-start+1)
			}
		} else {
			if err == io.EOF {
				break
			} else {
				log.Panic(err.Error())
			}
		}
	}
	return mrna
}

//func CreateMrnaFile(refgenePaths []string, genomeFastaPaths []string, mrnaFastaPath string) {
//	chromSeq := make(seq.Fasta)
//	for _, inGenomeFastaFile := range genomeFastaPaths {
//		for chrom, sequence := range seq.NewFasta(inGenomeFastaFile) {
//			chromSeq[chrom] = sequence
//		}
//	}
//	for _, refgenePath := range refgenePaths {
//		mrna := NewMrna(refgenePath, chromSeq)
//		mrna.Write(mrnaFastaPath)
//	}
//}
