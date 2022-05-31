package io

import (
	"fmt"
	"io"
	"open-anno/pkg/seq"
	"sort"
	"strings"

	"github.com/brentp/vcfgo"
)

type VCF struct {
	Chrom  string  `json:"CHROM"`
	Pos    int     `json:"POS"`
	ID     string  `json:"ID"`
	Ref    string  `json:"REF"`
	Alt    string  `json:"ALT"`
	Qual   float32 `json:"QUAL"`
	Filter string  `json:"FILTER"`
	Info   string  `json:"INFO"`
}

func (this VCF) Variant() Variant {
	chrom, pos, ref, alt := this.Chrom, this.Pos, this.Ref, this.Alt
	if chrom == "M" {
		chrom = "MT"
	}
	start, ref, alt := pos, strings.ToUpper(ref), strings.ToUpper(alt)
	if len(ref) > 1 || len(alt) > 1 && ref != alt {
		if strings.HasPrefix(ref, alt) || strings.HasSuffix(ref, alt) {
			if strings.HasPrefix(ref, alt) {
				start += len(alt)
			}
			ref = strings.Replace(ref, alt, "", 1)
			alt = ""
		} else if strings.HasPrefix(alt, ref) || strings.HasSuffix(alt, ref) {
			if strings.HasPrefix(alt, ref) {
				start += len(ref) - 1
			} else {
				start += len(ref) - len(alt)
			}
			alt = strings.Replace(alt, ref, "", 1)
			ref = ""
		} else {
			refRev := seq.Reverse(ref)
			altRef := seq.Reverse(alt)
			var length int
			length = seq.DifferenceSimple(refRev, altRef) - 1
			ref = ref[0 : len(ref)-length]
			alt = ref[0 : len(alt)-length]
			length = seq.DifferenceSimple(ref, alt) - 1
			ref = ref[length:]
			alt = alt[length:]
			start += length
			if length > 0 && len(ref) == 0 {
				start--
			}
		}
	}
	var end int
	if len(ref) == 0 {
		end = start
		ref = "-"
	} else {
		end = start + len(ref) - 1
	}
	if len(alt) == 0 {
		alt = "-"
	}
	return Variant{Chrom: chrom, Start: start, End: end, Ref: ref, Alt: alt, Otherinfo: []string{this.Info}}
}

type VCFs []VCF

func (this VCFs) Len() int           { return len(this) }
func (this VCFs) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this VCFs) Less(i, j int) bool { return this[i].Pos < this[j].Pos }

type VCFScanner struct {
	reader *vcfgo.Reader
	buffer *vcfgo.Variant
}

func NewVCFScanner(reader io.Reader) (VCFScanner, error) {
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return VCFScanner{}, err
	}
	return VCFScanner{reader: vcfReader}, err
}

func (this *VCFScanner) Close() error {
	return this.reader.Close()
}

func (this *VCFScanner) Scan() bool {
	row := this.reader.Read()
	if row == nil {
		return false
	}
	this.buffer = row
	return true
}

func (this VCFScanner) Text() string {
	return this.buffer.String()
}

func (this VCFScanner) Row() (VCFs, error) {
	vcfs := make(VCFs, len(this.buffer.Alt()))
	sample := this.buffer.Samples[0]
	depth := sample.DP
	alleles, err := sample.AltDepths()
	if err != nil {
		return vcfs, err
	}
	for i, alt := range this.buffer.Alt() {
		info := fmt.Sprintf("DEPTH=%d;VAF=%f;", depth, float64(alleles[i+1])/float64(depth))
		vcf := VCF{
			Chrom:  this.buffer.Chromosome,
			Pos:    int(this.buffer.Pos),
			ID:     this.buffer.Id(),
			Ref:    this.buffer.Reference,
			Alt:    alt,
			Qual:   this.buffer.Quality,
			Filter: this.buffer.Filter,
			Info:   info,
		}
		vcfs = append(vcfs, vcf)
	}
	return vcfs, err
}

func ReadVCFs(infile string) (VCFs, error) {
	var vcfs VCFs
	reader, err := NewIoReader(infile)
	if err != nil {
		return vcfs, err
	}
	defer reader.Close()
	scanner, err := NewVCFScanner(reader)
	if err != nil {
		return vcfs, err
	}
	defer scanner.Close()
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return vcfs, err
		}
		vcfs = append(vcfs, row...)
	}
	sort.Sort(vcfs)
	return vcfs, err
}

func WriteVCFs(outfile string, vcfs ...VCFs) error {
	writer, err := NewIoWriter(outfile)
	if err != nil {
		return err
	}
	fmt.Fprint(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for _, rows := range vcfs {
		for _, row := range rows {
			fmt.Fprintf(writer, "%s\t%d\t%s\t%s\t%s\t%f\t%s\t%s\n",
				row.Chrom, row.Pos, row.ID, row.Ref, row.Alt, row.Qual, row.Filter, row.Info)
		}
	}
	return nil
}
