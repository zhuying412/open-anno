package pkg

import (
	"fmt"
	"regexp"
)

func Min[T int | float64](a T, b T) T {
	if a < b {
		return a
	}
	return b
}

func Max[T int | float64](a T, b T) T {
	if a < b {
		return b
	}
	return a
}

func Abs[T int | float64](a T) T {
	if a < 0 {
		return -a
	}
	return a
}

func FindArr[T string | int | float64](arr []T, k T) int {
	for i, v := range arr {
		if k == v {
			return i
		}
	}
	return -1
}

func CurBin[T int | int64](chrom string, start T, size int) string {
	return fmt.Sprintf("%s:%d", chrom, start-(start%T(size)))
}

func FormatChrom(chrom string) string {
	chrom = regexp.MustCompile(`^chr`).ReplaceAllString(chrom, "")
	if chrom == "M" {
		return "MT"
	}
	return chrom
}
