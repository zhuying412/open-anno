package variant

import (
	"open-anno/pkg"
	"strconv"
	"strings"
)

type FilterBased struct {
	Chrom     string   `json:"chrom"`
	Start     int      `json:"start"`
	End       int      `json:"end"`
	Ref       string   `json:"ref"`
	Alt       string   `json:"alt"`
	Otherinfo []string `json:"otherinfo"`
}

func (this FilterBased) GetBaseVar() (string, int, int, string, string) {
	return this.Chrom, this.Start, this.End, this.Ref, this.Alt
}

func ReadFilterBasedLine(line string) (FilterBased, error) {
	fields := strings.Split(line, "\t")
	variant := FilterBased{
		Chrom:     pkg.FormatChrom(fields[0]),
		Ref:       fields[3],
		Alt:       fields[4],
		Otherinfo: fields[5:],
	}
	var err error
	variant.Start, err = strconv.Atoi(fields[1])
	if err != nil {
		return variant, err
	}
	variant.End, err = strconv.Atoi(fields[2])
	if err != nil {
		return variant, err
	}
	return variant, err
}
