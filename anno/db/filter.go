package db

import (
	"fmt"
	"log"
	"open-anno/pkg/io"
	"os"
	"sort"
	"strings"
)

func AnnoFilterBased(variants io.Variants, dbfile string, headerWrited bool, writer *os.File) {
	sort.Sort(variants)
	reader, err := os.Open(dbfile)
	if err != nil {
		log.Fatal(err)
	}
	defer reader.Close()
	scanner := io.NewDBVarScanner(reader)
	if headerWrited {
		writer.WriteString(scanner.Header + "\n")
	}
	var dbvar io.Variant
	for i := 0; i < len(variants); {
		switch variants[i].Compare(dbvar) {
		case io.VCMP_LT:
			i++
		case io.VCMP_GT:
			if !scanner.Scan() {
				break
			}
			dbvar, err = scanner.Row()
			if err != nil {
				log.Fatal(err)
			}
		case io.VCMP_EQ:
			writer.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s\t%s\t%s\n",
				dbvar.Chrom, dbvar.Start, dbvar.End,
				dbvar.Ref, dbvar.Alt, strings.Join(dbvar.Otherinfo, "\t"),
			))
			i++
		}
	}
}
