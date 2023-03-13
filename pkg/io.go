package pkg

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
	"strings"
)

func NewIOReader(infile string) (io.ReadCloser, error) {
	fi, err := os.Open(infile)
	if (strings.HasSuffix(strings.ToLower(infile), ".gz") || strings.HasSuffix(strings.ToLower(infile), ".bgz")) && err == nil {
		return gzip.NewReader(fi)
	}
	return fi, err
}

func NewIOWriter(infile string) (io.WriteCloser, error) {
	if strings.HasPrefix(infile, "-") {
		return os.Stdout, nil
	}
	fi, err := os.Create(infile)
	if strings.HasSuffix(strings.ToLower(infile), ".gz") && err == nil {
		return gzip.NewWriter(fi), nil
	}
	return fi, err
}

func IOCopy(src string, dst string) error {
	reader, err := os.Open(src)
	if err != nil {
		return err
	}
	defer reader.Close()
	writer, err := os.Create(dst)
	if err != nil {
		return err
	}
	if _, err := io.Copy(writer, reader); err != nil {
		return err
	}
	defer writer.Close()
	return err
}

type IOScanner struct {
	bufio.Scanner
}

func NewIOScanner(reader io.ReadCloser) IOScanner {
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 64*1024*1024), 1024*1024*1024)
	return IOScanner{Scanner: *scanner}
}

type CSVScanner struct {
	IOScanner
	FieldNames []string
}

func NewCSVScanner(reader io.ReadCloser) CSVScanner {
	scanner := NewIOScanner(reader)
	scanner.Scan()
	return CSVScanner{IOScanner: scanner, FieldNames: strings.Split(scanner.Text(), "\t")}
}

func (this CSVScanner) Row() map[string]string {
	row := make(map[string]string)
	fields := strings.Split(this.Text(), "\t")
	for i, key := range this.FieldNames {
		row[key] = fields[i]
	}
	return row
}
