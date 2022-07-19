package io

import (
	"bufio"
	"compress/gzip"
	"io"
	"log"
	"os"
	"strings"
)

type Reader io.Reader
type Writer io.Writer
type ReadCloser io.ReadCloser
type WriteCloser io.WriteCloser

const (
	SeekStart   = io.SeekStart
	SeekCurrent = io.SeekCurrent
	SeekEnd     = io.SeekEnd
)

var EOF = io.EOF

func NewIoReader(infile string) (ReadCloser, error) {
	fi, err := os.Open(infile)
	if strings.HasSuffix(strings.ToLower(infile), ".gz") && err == nil {
		return gzip.NewReader(fi)
	}
	return fi, err
}

func NewIoWriter(infile string) (WriteCloser, error) {
	if strings.HasPrefix(infile, "-") {
		return os.Stdout, nil
	}
	fi, err := os.Create(infile)
	if strings.HasSuffix(strings.ToLower(infile), ".gz") && err == nil {
		return gzip.NewWriter(fi), nil
	}
	return fi, err
}

func CopyFile(src string, dst string) error {
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

type IScanner[T any] interface {
	Scan() bool
	Text() string
	Row() (T, error)
}

type Scanner[T any] struct {
	bufio.Scanner
}

func (this Scanner[T]) Row() (T, error) {
	log.Fatal("Not implemented")
	return *new(T), nil
}

func NewScanner[T any](reader io.Reader) Scanner[T] {
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	return Scanner[T]{Scanner: *scanner}
}
