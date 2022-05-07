package database

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"os"
	"sort"
	"strings"
)

func AnnoRegionBased(variants io.Variants, dbfile string, overlap float64, writeHeader bool, writer *os.File) {
	sort.Sort(variants)
	regionBaseds, header, err := io.ReadDBBEDs(dbfile)
	if err != nil {
		log.Fatal(err)
	}
	if writeHeader {
		writer.WriteString(header + "\n")
	}
	for _, variant := range variants {
		var annos []string
		for _, dbvar := range regionBaseds {
			if variant.End >= dbvar.Start && variant.Start <= dbvar.End {
				vlen := variant.End - variant.Start + 1
				olen := pkg.Min(variant.End, dbvar.End) - pkg.Max(variant.Start, dbvar.Start) + 1
				if float64(olen)/float64(vlen) >= overlap {
					annos = append(annos, dbvar.Name)
				}
			}
		}
		if len(annos) > 0 {
			writer.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s\t%s\t%s\n",
				variant.Chrom, variant.Start, variant.End,
				variant.Ref, variant.Alt, strings.Join(annos, ","),
			))
		}
	}
}
