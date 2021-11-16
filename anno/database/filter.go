package database

import (
	"OpenAnno/variant"
	"bufio"
	"log"
	"os"
)

func RunFilterAnnotate(variants variant.IVariants, databaseFile string) map[string]DBAnno {
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
	var _variant variant.Variant
	var _anno DBAnno
	var headers []string
	annoMap := make(map[string]DBAnno)
	for i := 0; i < variants.Len(); {
		if ReadLine(reader, &headers, &_variant, &_anno) != nil {
			break
		}
		cmp := variants.GetVariant(i).Compare(_variant.Start, _variant.End, _variant.Ref, _variant.Alt)
		if cmp == variant.VarCmp_LT {
			i++
		} else if cmp == variant.VarCmp_GT {
			if ReadLine(reader, &headers, &_variant, &_anno) != nil {
				break
			}
		} else {
			if cmp == variant.VarCmp_EQ {
				annoMap[variants.GetVariant(i).SN()] = _anno
				i++
			} else {
				if ReadLine(reader, &headers, &_variant, &_anno) != nil {
					break
				}
			}
		}
	}
	return annoMap
}
