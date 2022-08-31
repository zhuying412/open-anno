package db

import (
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"open-anno/pkg/scheme"
	"os"
	"sort"
)

func ReadFilterBasedIndexMap(infile string) (map[string]scheme.DBVarIdx, int, error) {
	idxMap := make(map[string]scheme.DBVarIdx)
	reader, err := io.NewIoReader(infile)
	if err != nil {
		return idxMap, 0, err
	}
	scanner := io.NewDBVarIdxScanner(reader)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return idxMap, scanner.BinSize, err
		}
		idxMap[row.Bin] = row
	}
	return idxMap, scanner.BinSize, err
}

func ReadVariantMap(avinput string, binSize int) (map[string]scheme.Variants, error) {
	variants := make(map[string]scheme.Variants)
	reader, err := io.NewIoReader(avinput)
	if err != nil {
		return variants, err
	}
	defer reader.Close()
	scanner := io.NewVarScanner(reader)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return variants, err
		}
		curbin := pkg.CurBin(row.Chrom, row.Start, binSize)
		if rows, ok := variants[curbin]; ok {
			variants[curbin] = append(rows, row)
		} else {
			variants[curbin] = scheme.Variants{row}
		}
	}
	return variants, err
}

func AnnoFilterBased(infile, dbfile, outfile string) error {
	// 读取 index file
	idxMap := make(map[string]scheme.DBVarIdx)
	idxReader, err := io.NewIoReader(dbfile + ".idx")
	if err != nil {
		return err
	}
	defer idxReader.Close()
	idxScanner := io.NewDBVarIdxScanner(idxReader)
	for idxScanner.Scan() {
		row, err := idxScanner.Row()
		if err != nil {
			return err
		}
		idxMap[row.Bin] = row
	}
	binSize := idxScanner.BinSize
	// 读取变异
	variantBinMap, err := ReadVariantMap(infile, binSize)
	if err != nil {
		return err
	}
	// 定义输出
	writer, err := io.NewIoWriter(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	// 开始注释
	reader, err := os.Open(dbfile)
	if err != nil {
		return err
	}
	defer reader.Close()
	scanner := io.NewDBVarScanner(reader)
	fmt.Fprintf(writer, "%s\n", scanner.Header)
	for curbin, variants := range variantBinMap {
		sort.Sort(variants)
		if idx, ok := idxMap[curbin]; ok {
			reader.Seek(idx.Start, io.SeekStart)
			var dbvar scheme.DBVar
			for i, length := 0, idx.Start; i < len(variants) && length <= idx.End; {
				switch dbvar.Compare(variants[i]) {
				case scheme.VCMP_GT:
					i++
				case scheme.VCMP_LT:
					if !scanner.Scan() {
						break
					}
					length += int64(len(scanner.Text()))
					dbvar, err = scanner.Row()
					if err != nil {
						return err
					}
				case scheme.VCMP_EQ:
					fmt.Fprintf(writer, "%s\n", dbvar.Text)
					i++
				}
			}

		}
	}
	return err
}
