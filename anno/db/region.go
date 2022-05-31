package db

import (
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"sort"
	"strings"
)

func AnnoRegionBased(variants io.Variants, dbfile string, overlap float64, writeHeader bool, writer io.WriteCloser) error {
	sort.Sort(variants)
	regionBaseds, header, err := io.ReadDBBEDs(dbfile)
	if err != nil {
		return err
	}
	if writeHeader {
		headers := strings.Split(header, "\t")
		fmt.Fprintf(writer, "%s\t%s\t%s\tRef\tAlt\t%s\n", headers[0], headers[1], headers[2], headers[3])
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
			fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\n",
				variant.Chrom, variant.Start, variant.End,
				variant.Ref, variant.Alt, strings.Join(annos, ","))
		}
	}
	return err
}
