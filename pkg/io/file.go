package io

import (
	"compress/gzip"
	"io"
	"os"
	"strings"
)

type Reader io.Reader
type Writer io.Writer
type ReadCloser io.ReadCloser
type WriteCloser io.WriteCloser

func NewIoReader(infile string) (ReadCloser, error) {
	fi, err := os.Open(infile)
	if strings.HasSuffix(strings.ToLower(infile), ".gz") && err == nil {
		return gzip.NewReader(fi)
	}
	return fi, err
}

func NewIoWriter(infile string) (WriteCloser, error) {
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
