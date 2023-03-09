package anno

import (
	"errors"
	"fmt"
	"open-anno/pkg"
	"strconv"
	"strings"

	"github.com/brentp/vcfgo"
)

const (
	VType_SNP = "SNP"
	VType_INS = "INS"
	VType_DEL = "DEL"
	VType_DUP = "DUP"
	VType_SUB = "SUB"
)

type AnnoVariant struct {
	Chrom string `json:"chrom"`
	Start int    `json:"start"`
	End   int    `json:"end"`
	Ref   string `json:"ref"`
	Alt   string `json:"alt"`
}

func (this AnnoVariant) Type() string {
	if this.Ref == "DIP" {
		return this.Alt
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

func (this AnnoVariant) PK() string {
	return fmt.Sprintf("%s:%d:%d:%s:%s", this.Chrom, this.Start, this.End, this.Ref, this.Alt)
}

func (this AnnoVariant) Less(that AnnoVariant) bool {
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

func (this AnnoVariant) Equal(that AnnoVariant) bool {
	return this.Chrom == that.Chrom && this.Start == that.Start && this.End == that.End && this.Ref == that.Ref && this.Alt == that.Alt
}

func (this AnnoVariant) Greater(that AnnoVariant) bool {
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

type IVariant interface {
	AnnoVariant() AnnoVariant
}

type Variants []IVariant

func (this Variants) Len() int           { return len(this) }
func (this Variants) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this Variants) Less(i, j int) bool { return this[i].AnnoVariant().Less(this[j].AnnoVariant()) }

func (this Variants) AggregateByChrom() map[string]Variants {
	variants := make(map[string]Variants)
	for _, variant := range this {
		chrom := variant.AnnoVariant().Chrom
		if rows, ok := variants[chrom]; ok {
			variants[chrom] = append(rows, variant)
		} else {
			variants[chrom] = Variants{variant}
		}
	}
	return variants
}

func (this Variants) AggregateByBin(binSize int) map[string]Variants {
	variants := make(map[string]Variants)
	for _, variant := range this {
		chrom, start := variant.AnnoVariant().Chrom, variant.AnnoVariant().Start
		curbin := pkg.CurBin(chrom, start, binSize)
		if rows, ok := variants[curbin]; ok {
			variants[curbin] = append(rows, variant)
		} else {
			variants[curbin] = Variants{variant}
		}
	}
	return variants
}

type SNV struct {
	vcfgo.Variant
}

func (this SNV) AnnoVariant() AnnoVariant {
	alt := this.Alt()[0]
	chrom, start, end, ref, alt := pkg.VCFtoAV(this.Chrom(), int(this.Pos), this.Ref(), alt)
	return AnnoVariant{Chrom: chrom, Start: start, End: end, Ref: ref, Alt: alt}
}

func ReadVCF(vcfFile string) (Variants, *vcfgo.Header, error) {
	snvs := make(Variants, 0)
	reader, err := pkg.NewIOReader(vcfFile)
	if err != nil {
		return snvs, &vcfgo.Header{}, err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return snvs, &vcfgo.Header{}, err
	}
	defer vcfReader.Close()
	for variant := vcfReader.Read(); variant != nil; variant = vcfReader.Read() {
		if len(variant.Alt()) > 1 {
			return snvs, &vcfgo.Header{}, errors.New(fmt.Sprintf("the count of alt > 1: %s", variant.String()))
		}
		snvs = append(snvs, SNV{Variant: *variant})
	}
	return snvs, vcfReader.Header, nil
}

type CNV struct {
	Chrom string `json:"chrom"`
	Start int    `json:"start"`
	End   int    `json:"end"`
	Type  string `json:"type"`
}

func (this CNV) AnnoVariant() AnnoVariant {
	alt := strings.ToUpper(this.Type)
	if alt == "DUP" || alt == "DUPLICATION" || alt == "GAIN" {
		alt = "DUP"
	} else if alt == "DEL" || alt == "DELETION" || alt == "LOSS" {
		alt = "DEL"
	}
	return AnnoVariant{Chrom: this.Chrom, Start: this.Start, End: this.End, Ref: "DIP", Alt: alt}
}

func ReadBED(bedFile string) (Variants, error) {
	cnvs := make(Variants, 0)
	reader, err := pkg.NewIOReader(bedFile)
	if err != nil {
		return cnvs, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		start, err := strconv.Atoi(row[1])
		if err != nil {
			return cnvs, err
		}
		end, err := strconv.Atoi(row[2])
		if err != nil {
			return cnvs, err
		}
		cnv := CNV{
			Chrom: row[0],
			Start: start,
			End:   end,
			Type:  row[3],
		}
		cnvs = append(cnvs, cnv)
	}
	return cnvs, nil
}
