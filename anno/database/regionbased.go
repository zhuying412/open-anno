package database

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/variant"
	"os"
	"sort"
	"strings"
)

func AnnoRegionBased(variants variant.Variants, dbfile string, overlap float64, writer *os.File) {
	sort.Sort(variants)
	regionBaseds, header, err := variant.ReadRegionBasedDB(dbfile)
	if err != nil {
		log.Fatal(err)
	}
	writer.WriteString(header + "\n")
	for _, variant := range variants {
		var info []string
		for _, dbvar := range regionBaseds {
			if variant.End >= dbvar.Start && variant.Start <= dbvar.End {
				vlen := variant.End - variant.Start + 1
				olen := pkg.Min(variant.End, dbvar.End) - pkg.Max(variant.Start, dbvar.Start) + 1
				if float64(olen)/float64(vlen) >= overlap {
					info = append(info, dbvar.Otherinfo)
				}
			}
		}
		anno := strings.Join(info, ",")
		if anno == "" {
			anno = "."
		}
		writer.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s\t%s\t%s\n",
			variant.Chrom, variant.Start, variant.End,
			variant.Ref, variant.Alt, anno,
		))
	}
}
