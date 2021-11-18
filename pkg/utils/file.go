package utils

import (
	"bufio"
	"io"
	"log"
	"os"
	"strings"
)

func OpenFile(filePath string) (*os.File, *bufio.Reader) {
	fi, err := os.Open(filePath)
	if err != nil {
		log.Panicln(fi)
	}
	reader := bufio.NewReader(fi)
	return fi, reader
}

func CloseFile(fi *os.File) {
	err := fi.Close()
	if err != nil {
		log.Panic(err)
	}
}

func ReadLine(reader *bufio.Reader, comment byte) (string, bool) {
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				return "", true
			}
		}
		line = strings.TrimSpace(line)
		if len(line) == 0 {
			continue
		}
		if comment != 0 && line[0] == comment {
			continue
		}
		return line, false
	}
}

func WriteLine(fo *os.File, line string) {
	_, err := fo.WriteString(line)
	if err != nil {
		log.Panic(err)
	}
}

func CreateFile(filePath string) *os.File {
	fo, err := os.Create(filePath)
	if err != nil {
		log.Panic(err)
	}
	return fo
}

func CreateDir(dirPath string) {
	err := os.MkdirAll(dirPath, os.ModePerm)
	if err != nil {
		log.Panic(err)
	}
}
