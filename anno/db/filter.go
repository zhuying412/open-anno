package db

import (
	"fmt"
	"log"
	"open-anno/pkg/io"
	"os"
	"sort"
	"strings"
)

func AnnoFilterBased(variants io.Variants, dbfile string, headerWrited bool, writer io.WriteCloser) {
	sort.Sort(variants)
	reader, err := os.Open(dbfile)
	if err != nil {
		log.Fatal(err)
	}
	defer reader.Close()
	scanner := io.NewDBVarScanner(reader)
	if headerWrited {
		fmt.Fprintf(writer, "%s\n", scanner.Header)
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
			fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\n",
				dbvar.Chrom, dbvar.Start, dbvar.End,
				dbvar.Ref, dbvar.Alt, strings.Join(dbvar.Otherinfo, "\t"))
			i++
		}
	}
}
