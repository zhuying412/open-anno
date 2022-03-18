package variant

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

const (
	VType_SNP = "SNP"
	VType_INS = "INS"
	VType_DEL = "DEL"
	VType_DUP = "DUP"
)

type Variant struct {
	Chrom     string   `json:"chrom"`
	Start     int      `json:"start"`
	End       int      `json:"end"`
	Ref       string   `json:"ref"`
	Alt       string   `json:"alt"`
	Otherinfo []string `json:"otherinfo"`
}

func (this Variant) Type() string {
	if this.Ref == "DIP" {
		if this.Alt == "DEL" {
			return VType_DEL
		}
		return VType_DUP
	}
	if this.Alt == "-" {
		return VType_DEL
	} else if this.Ref == "-" {
		return VType_INS
	}
	return VType_SNP
}

func (this Variant) ID() string {
	return fmt.Sprintf("%s:%d-%d:%s/%s", this.Chrom, this.Start, this.End, this.Ref, this.Alt)
}

type Variants []Variant

func (this Variants) Len() int           { return len(this) }
func (this Variants) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this Variants) Less(i, j int) bool { return this[i].Start < this[j].Start }

func ReadAvinput(avinput string) (map[string]Variants, error) {
	variants := make(map[string]Variants)
	fi, err := os.Open(avinput)
	if err != nil {
		return variants, err
	}
	defer fi.Close()
	scanner := bufio.NewScanner(fi)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		chrom := fields[0]
		if chrom[0] == '#' {
			continue
		}
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return variants, err
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return variants, err
		}
		variant := Variant{
			Chrom:     fields[0],
			Start:     start,
			End:       end,
			Ref:       fields[3],
			Alt:       fields[4],
			Otherinfo: fields[5:],
		}
		if vars, ok := variants[chrom]; ok {
			variants[chrom] = append(vars, variant)
		} else {
			variants[chrom] = Variants{variant}
		}
	}
	return variants, err
}
