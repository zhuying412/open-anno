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
