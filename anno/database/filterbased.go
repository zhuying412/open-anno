package database

import (
	"bufio"
	"fmt"
	"log"
	"open-anno/pkg/variant"
	"os"
	"sort"
	"strings"
)

func AnnoFilterBased(variants variant.Variants, dbfile string, headerWrited bool, writer *os.File) {
	sort.Sort(variants)
	fi, err := os.Open(dbfile)
	if err != nil {
		log.Fatal(err)
	}
	defer fi.Close()
	scanner := bufio.NewScanner(fi)
	scanner.Scan()
	if headerWrited {
		writer.WriteString(scanner.Text() + "\n")
	}
	var dbvar variant.FilterBased
	for i := 0; i < len(variants); {
		switch variant.CompareVar(variants[i], dbvar) {
		case variant.VCMP_LT:
			i++
		case variant.VCMP_GT:
			if !scanner.Scan() {
				break
			}
			dbvar, err = variant.ReadFilterBasedLine(scanner.Text())
			if err != nil {
				log.Fatal(err)
			}
		case variant.VCMP_EQ:
			writer.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s\t%s\t%s\n",
				dbvar.Chrom, dbvar.Start, dbvar.End,
				dbvar.Ref, dbvar.Alt, strings.Join(dbvar.Otherinfo, "\t"),
			))
			i++
		}
	}
}
