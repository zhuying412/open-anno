package tools

import (
	"bufio"
	"fmt"
	"open-anno/pkg/variant"
	"os"
	"strings"
)

type BED struct {
	Chrom string `json:"chrom"`
	Start int    `json:"start"`
	END   int    `json:"end"`
	Type  string `json:"type"`
}

type BEDs []BED

func (this BEDs) Len() int           { return len(this) }
func (this BEDs) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this BEDs) Less(i, j int) bool { return this[i].Start < this[j].Start }

func ReadCnvAV(infile string) (BEDs, error) {
	var beds BEDs
	fi, err := os.Open(infile)
	if err != nil {
		return beds, err
	}
	defer fi.Close()
	scanner := bufio.NewScanner(fi)
	for scanner.Scan() {
		line := scanner.Text()
		if line[0] == '#' {
			continue
		}
		variant, err := variant.ReadVariantLine(line)
		if err != nil {
			return beds, err
		}
		chrom, start, end := variant.Chrom, variant.Start, variant.End
		infoMap := make(map[string]string)
		for _, info := range variant.Otherinfo {
			for _, item := range strings.Split(info, ";") {
				if strings.Contains(item, "=") {
					kv := strings.Split(item, "=")
					infoMap[kv[0]] = kv[1]
				}
			}
		}
		if alt, ok := infoMap["ALT"]; ok {
			for _, typo := range strings.Split(alt, "/") {
				beds = append(beds, BED{Chrom: chrom, Start: start, END: end, Type: typo})
			}
		}

	}
	return beds, err
}

func WriteBED(beds BEDs, outfile string) error {
	writer, err := os.Create(outfile)
	if err != nil {
		return err
	}
	for _, bed := range beds {
		fmt.Fprintf(writer, "%s\t%d\t%d\t%s\n", bed.Chrom, bed.Start, bed.END, bed.Type)
	}
	return nil
}
