package anno

import (
	"fmt"
	"open-anno/pkg"
	"strconv"
	"strings"
)

const (
	VType_SNP  = "SNP"
	VType_INS  = "INS"
	VType_DEL  = "DEL"
	VType_DUP  = "DUP"
	VType_SUB  = "SUB"
	VType_GAIN = "GAIN"
	VType_LOSS = "LOSS"
)

type Variant struct {
	Chrom     string `json:"chrom"`
	Start     int    `json:"start"`
	End       int    `json:"end"`
	Ref       string `json:"ref"`
	Alt       string `json:"alt"`
	Otherinfo string `json:"otherinfo"`
}

func (this Variant) Type() string {
	if this.Ref == "DIP" {
		if this.Alt == "DEL" || this.Alt == "GAIN" {
			return VType_GAIN
		}
		return VType_LOSS
	}
	if this.Ref == "-" {
		return VType_INS
	} else if this.Alt == "-" {
		return VType_DEL
	} else {
		if len(this.Ref) > 1 || len(this.Alt) > 1 {
			return VType_SUB
		}
		return VType_SNP
	}
}

func (this Variant) PK() string {
	return fmt.Sprintf("%s:%d:%d:%s:%s", this.Chrom, this.Start, this.End, this.Ref, this.Alt)
}

func (this Variant) Less(that Variant) bool {
	if this.Chrom < that.Chrom {
		return true
	}
	if this.Chrom == that.Chrom {
		if this.Start < that.Start {
			return true
		}
		if this.Start == that.Start {
			if this.End < that.End {
				return true
			}
			if this.End == that.End {
				if this.Ref < that.Ref {
					return true
				}
				if this.Ref == that.Ref {
					if this.Alt < that.Alt {
						return true
					}
				}
			}
		}
	}
	return false
}

func (this Variant) Equal(that Variant) bool {
	return this.Chrom == that.Chrom && this.Start == that.Start && this.End == that.End && this.Ref == that.Ref && this.Alt == that.Alt
}

func (this Variant) Greater(that Variant) bool {
	if this.Chrom > that.Chrom {
		return true
	}
	if this.Chrom == that.Chrom {
		if this.Start > that.Start {
			return true
		}
		if this.Start == that.Start {
			if this.End > that.End {
				return true
			}
			if this.End == that.End {
				if this.Ref > that.Ref {
					return true
				}
				if this.Ref == that.Ref {
					if this.Alt > that.Alt {
						return true
					}
				}
			}
		}
	}
	return false
}

type Variants []Variant

func (this Variants) Len() int           { return len(this) }
func (this Variants) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this Variants) Less(i, j int) bool { return this[i].Less(this[j]) }

func ReadAnnoInput(annoInputFile string) (map[string]Variants, error) {
	variants := make(map[string]Variants)
	reader, err := pkg.NewIOReader(annoInputFile)
	if err != nil {
		return variants, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		variant := Variant{
			Chrom:     row[0],
			Ref:       row[3],
			Alt:       row[4],
			Otherinfo: row[5],
		}
		variant.Start, err = strconv.Atoi(row[1])
		if err != nil {
			return variants, err
		}
		variant.End, err = strconv.Atoi(row[2])
		if err != nil {
			return variants, err
		}
		if rows, ok := variants[variant.Chrom]; ok {
			variants[variant.Chrom] = append(rows, variant)
		} else {
			variants[variant.Chrom] = Variants{variant}
		}
	}
	return variants, err
}
