package pkg

import (
	"fmt"
	"os"

	"github.com/go-playground/validator/v10"
)

func Min[T int | float64 | uint64 | uint32](a T, b T) T {
	if a < b {
		return a
	}
	return b
}

func Max[T int | float64 | uint64 | uint32](a T, b T) T {
	if a < b {
		return b
	}
	return a
}

func Abs[T int | float64 | uint64 | uint32](a T) T {
	if a < 0 {
		return -a
	}
	return a
}

func Sum[T int | float64 | uint64 | uint32](a ...T) T {
	var sum T
	for _, i := range a {
		sum += i
	}
	return sum
}

func FindArr[T string | int | float64 | uint64 | uint32](arr []T, k T) int {
	for i, v := range arr {
		if k == v {
			return i
		}
	}
	return -1
}

func Interface2Array[T string | int | float64 | uint64 | uint32](src any) []T {
	if dst, ok := src.(T); ok {
		return []T{dst}
	}
	return src.([]T)
}

func UniqArr[T string | int | float64 | uint64 | uint32](arr []T, ex T) []T {
	uarr := make([]T, 0)
	for _, v := range arr {
		if v != ex && FindArr(arr, v) < 0 {
			uarr = append(uarr, v)
		}
	}
	return uarr
}

func SplitArr[T any](arr []T, num int) [][]T {
	length := len(arr)
	size := length / num
	if len(arr)%num != 0 {
		size++
	}
	var segments = make([][]T, 0)
	var start, end int
	for i := 1; i <= num; i++ {
		end = i * size
		if end > length {
			segments = append(segments, arr[start:])
			break
		}
		segments = append(segments, arr[start:end])
		start = i * size
	}
	return segments
}

func NewArr[T any](size int, fill T) []T {
	arr := make([]T, size)
	for i := 0; i < size; i++ {
		arr[i] = fill
	}
	return arr
}

func CheckPathExists(fl validator.FieldLevel) bool {
	path := fl.Field().String()
	_, err := os.Stat(path)
	return !os.IsNotExist(err)
}

func CheckPathsExists(fl validator.FieldLevel) bool {
	for _, path := range fl.Field().Interface().([]string) {
		_, err := os.Stat(path)
		if os.IsNotExist(err) {
			return false
		}
	}
	return true
}

// CurBin 计算Bin游标
func CurBin[T int | int64](chrom string, start T, size int) string {
	return fmt.Sprintf("%s\t%d", chrom, start-(start%T(size)))
}
