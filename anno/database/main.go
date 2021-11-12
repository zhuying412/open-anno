package database

import (
	"OpenAnno/variant"
	"bufio"
	"log"
	"os"
)

func RunAnnotate(variants variant.IVariants, databaseFile string) map[string]IAnno {
	fi, err := os.Open(databaseFile)
	if err != nil {
		log.Panic(err)
	}
	defer func(fp *os.File) {
		err = fp.Close()
		if err != nil {
			log.Panic(err)
		}
	}(fi)
	reader := bufio.NewReader(fi)
	annoMap := make(map[string]IAnno)
	var headers []string
	var _variant variant.Variant
	var _anno IAnno
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
			_anno.ProcessOverlap()
		}
	}
	return annoMap
}
