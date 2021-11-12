package region

import (
	"OpenAnno/variant"
	"bufio"
	"io"
	"log"
	"os"
	"strings"
)

func ReadLine(reader *bufio.Reader) (line string, isEof bool) {
	line, err := reader.ReadString('\n')
	if err != nil {
		if err == io.EOF {
			return line, true
		} else {
			log.Panic(err)
		}
	}
	return line, false
}

func RunAnnotate(variants variant.IVariants, databaseFile string) map[string]RegionBasedAnno {
	fi, err := os.Open(databaseFile)
	if err != nil {
		log.Panic(err)
	}
	defer func(fp *os.File) {
		err := fp.Close()
		if err != nil {
			log.Panic(err)
		}
	}(fi)
	reader := bufio.NewReader(fi)
	annoMap := make(map[string]RegionBasedAnno)
	var headers []string
	var _variant variant.Variant
	var _anno RegionBasedAnno
	for i := 0; i < variants.Len(); {
		if len(headers) == 0 {
			if header, isEof := ReadLine(reader); isEof {
				break
			} else {
				if len(header) == 0 || header[0] == '#' {
					continue
				}
				headers = strings.Split(header, "\t")
			}
		}
		if len(_anno) == 0 {
			if line, isEof := ReadLine(reader); isEof {
				break
			} else {
				if len(line) == 0 || line[0] == '#' {
					continue
				}
				_variant, _anno = ReadRegionFields(headers, strings.Split(line, "\t"))
			}
		}
		cmp := variants.GetVariant(i).Compare(_variant.Start, _variant.End, _variant.Ref, _variant.Alt)
		if cmp == variant.VarCmp_LT {
			i++
		} else if cmp == variant.VarCmp_GT {
			if line, isEof := ReadLine(reader); isEof {
				break
			} else {
				if len(line) == 0 || line[0] == '#' {
					continue
				}
				_variant, _anno = ReadRegionFields(headers, strings.Split(line, "\t"))
			}
		} else {
			annoMap[variants.GetVariant(i).SN()] = _anno
			i++
		}
	}
	return annoMap
}
