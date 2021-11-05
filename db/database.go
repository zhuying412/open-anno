package db

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

type DatabaseIndex struct {
	Chrom      string
	ChromStart int
	ChromEnd   int
	LineStart  int
	LineEnd    int
}

type DatabaseIndexes []DatabaseIndex

func (f DatabaseIndexes) FilterByChrom(chrom string) DatabaseIndexes {
	indexes := make(DatabaseIndexes, 0)
	for _, index := range f {
		if index.Chrom == chrom {
			indexes = append(indexes, index)
		}
	}
	return indexes
}

func ReadDatabaseIndexFile(DatabaseIndexFile string) DatabaseIndexes {
	indexes := make(DatabaseIndexes, 0)
	if fp, err := os.Open(DatabaseIndexFile); err == nil {
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
				chromStart, err := strconv.Atoi(fields[1])
				if err != nil {
					log.Panic(err)
				}
				chromEnd, err := strconv.Atoi(fields[2])
				if err != nil {
					log.Panic(err)
				}
				lineStart, err := strconv.Atoi(fields[3])
				if err != nil {
					log.Panic(err)
				}
				lineEnd, err := strconv.Atoi(fields[4])
				if err != nil {
					log.Panic(err)
				}
				indexes = append(indexes, DatabaseIndex{
					Chrom:      fields[0],
					ChromStart: chromStart,
					ChromEnd:   chromEnd,
					LineStart:  lineStart,
					LineEnd:    lineEnd,
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
	return indexes
}

func ReadDatabaseFile(DatabaseFile string) DatabaseIndexes {
	indexes := make(DatabaseIndexes, 0)
	if fp, err := os.Open(DatabaseFile); err == nil {
		defer func(fp *os.File) {
			err := fp.Close()
			if err != nil {
				log.Panic(err)
			}
		}(fp)
		reader := bufio.NewReader(fp)
		lineCount, chromStart, chromEnd, chrom := 0, 0, 0, ""
		for {
			if line, err := reader.ReadString('\n'); err == nil {
				line = strings.TrimSpace(line)
				if len(line) > 0 && line[0] != '#' {
					fields := strings.Split(line, "\t")
					if (chrom != "" && fields[0] != chrom) ||
						(lineCount > 0 && lineCount%DBIndexStepLen == 0) {
						indexes = append(indexes, DatabaseIndex{
							Chrom:      chrom,
							ChromStart: chromStart,
							ChromEnd:   chromEnd,
							LineStart:  lineCount - DBIndexStepLen + 1,
							LineEnd:    lineCount,
						})
						//log.Printf("%s\t%d\t%d\t%d\t%d\n", chrom, chromStart, chromEnd, lineCount-DBIndexStepLen+1, lineCount)
						chromStart, chromEnd = 0, 0
						if fields[0] != chrom {
							lineCount = 0
						}
					}
					chrom = fields[0]
					start, err := strconv.Atoi(fields[1])
					if err != nil {
						log.Panic(err)
					}
					end, err := strconv.Atoi(fields[2])
					if err != nil {
						log.Panic(err)
					}
					if len(indexes) > 0 &&
						indexes[len(indexes)-1].Chrom == chrom &&
						indexes[len(indexes)-1].ChromEnd >= end {
						indexes[len(indexes)-1].LineEnd++
					} else {
						if chromStart == 0 {
							chromStart = start
						}
						chromEnd = end
					}
				}
				lineCount++
			} else {
				if err == io.EOF {
					break
				} else {
					log.Panic(err)
				}
			}
		}
		lineStart := lineCount - DBIndexStepLen + 1
		if lineStart < 0 {
			lineStart = 1
		}
		indexes = append(indexes, DatabaseIndex{
			Chrom:      chrom,
			ChromStart: chromStart,
			ChromEnd:   chromEnd,
			LineStart:  lineStart,
			LineEnd:    lineCount,
		})
		//log.Printf("%s\t%d\t%d\t%d\t%d\n", chrom, chromStart, chromEnd, lineStart, lineCount)
	} else {
		log.Panic(err)
	}
	return indexes
}

func CreateDatabseIndexFile(indexes DatabaseIndexes, DatabaseIndexFile string) {
	fo, err := os.Create(DatabaseIndexFile)
	if err == nil {
		defer func(fo *os.File) {
			err := fo.Close()
			if err != nil {
				log.Panic(err)
			}
		}(fo)
		for _, index := range indexes {
			if _, err := fo.WriteString(fmt.Sprintf(
				"%s\t%d\t%d\t%d\t%d\n",
				index.Chrom, index.ChromStart, index.ChromEnd,
				index.LineStart, index.LineEnd,
			)); err != nil {
				log.Panic(err)
			}
		}
	}
}
