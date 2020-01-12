package snv

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"grandanno/core"
	"io/ioutil"
	"os"
	"sort"
	"strconv"
	"strings"
)

type Information struct {
	Depth      int
	Qual       float32
	GatkFilter string
	Genotype   float32
	Ratio      float32
}

type Snv struct {
	Variant     core.Variant
	Information Information
	OtherInfo   string
}

type Snvs []Snv

func (snv Snv) GetType() string {
	var typo string
	switch {
	case snv.Variant.Ref.GetChar(0) == '-':
		typo = "ins"
	case snv.Variant.Alt.GetChar(0) == '-':
		typo = "del"
	default:
		typo = "snp"
	}
	return typo
}

func (snvs Snvs) Len() int {
	return len(snvs)
}

func (snvs Snvs) Less(i, j int) bool {
	starti, endi := snvs[i].Variant.GetDigitalPosition()
	startj, endj := snvs[i].Variant.GetDigitalPosition()
	if starti == startj {
		return endi < endj
	} else {
		return starti < startj
	}
}

func (snvs Snvs) Swap(i, j int) {
	snvs[i], snvs[j] = snvs[j], snvs[i]
}

type Vcf struct {
	File string
}

func (vcf Vcf) getGenotypes(info string, n int) []float32 {
	genotypes := make([]float32, n)
	for i := 0; i < n; i++ {
		genotypes[i] = -1
	}
	for _, s := range strings.Split(info, ";") {
		if strings.HasPrefix(s, "AF=") {
			for i, _gt := range strings.Split(s[3:], ",") {
				if i > n {
					break
				}
				if gt, err := strconv.ParseFloat(_gt, 32); err == nil {
					genotypes[i] = float32(gt)
				}
			}
			break
		}
	}
	return genotypes
}

func (vcf Vcf) getDepthRatios(formatKey string, formatValue string, n int) (int, []float32) {
	ratios := make([]float32, n)
	depth := -1
	for i := 0; i < n; i++ {
		ratios[i] = float32(-1)
	}
	keys := strings.Split(formatKey, ":")
	values := strings.Split(formatValue, ":")
	for i, key := range keys {
		if key == "DP" {
			if d, err := strconv.Atoi(values[i]); err == nil {
				depth = d
			}
		}
		if key == "AD" {
			sum := 0
			counts := make([]int, n)
			_counts := strings.Split(values[i], ",")
			for i := 0; i <= n && i < len(_counts); i++ {
				if c, err := strconv.Atoi(_counts[i]); err == nil {
					sum += c
					if i > 0 {
						counts[i-1] = c
					}
				}
			}
			if sum > 0 {
				for i := 0; i < n; i++ {
					ratios[i] = float32(counts[i]) / float32(sum)
				}
			}
		}
	}
	return depth, ratios
}

func (vcf Vcf) readLine(vcfLine string) Snvs {
	field := strings.Split(vcfLine, "\t")
	alts := strings.Split(field[4], ",")
	genotypes := vcf.getGenotypes(field[7], len(alts))
	depth, ratios := vcf.getDepthRatios(field[len(field)-2], field[len(field)-1], len(alts))
	gatkFilter := field[6]
	qual := float32(-1)
	if q, err := strconv.ParseFloat(field[5], 32); err == nil {
		qual = float32(q)
	}
	if pos, err := strconv.Atoi(field[1]); err == nil {
		snvs := make(Snvs, len(alts))
		for i, alt := range alts {
			snvs[i].Variant = core.Variant{
				Chrom: strings.Replace(field[0], "chr", "", 1),
				Start: pos,
				End:   0,
				Ref:   core.Sequence(field[3]),
				Alt:   core.Sequence(alt),
			}
			snvs[i].Variant.ConvertSnv()
			snvs[i].Information = Information{
				Depth:      depth,
				Qual:       qual,
				GatkFilter: gatkFilter,
				Genotype:   genotypes[i],
				Ratio:      ratios[i],
			}
			snvs[i].OtherInfo = vcfLine
		}
		return snvs
	}
	return Snvs{}
}

func (vcf Vcf) ReadAll() Snvs {
	var snvs Snvs
	if fp, err := os.Open(vcf.File); err == nil {
		defer fp.Close()
		var lines [][]byte
		if strings.HasSuffix(strings.ToLower(vcf.File), ".gz") {
			if reader, _err := gzip.NewReader(fp); _err == nil {
				if content, err := ioutil.ReadAll(reader); err == nil {
					lines = bytes.Split(content, []byte{'\n'})
				}
			}
		} else {
			reader := bufio.NewReader(fp)
			if content, err := ioutil.ReadAll(reader); err == nil {
				lines = bytes.Split(content, []byte{'\n'})
			}
		}
		for _, line := range lines {
			line = bytes.TrimSpace(line)
			if len(line) == 0 || line[0] == '#' {
				continue
			}
			field := strings.Split(string(line), "\t")
			alts := strings.Split(field[4], ",")
			genotypes := vcf.getGenotypes(field[7], len(alts))
			depth, ratios := vcf.getDepthRatios(field[len(field)-2], field[len(field)-1], len(alts))
			gatkFilter := field[6]
			qual := float32(-1)
			if q, err := strconv.ParseFloat(field[5], 32); err == nil {
				qual = float32(q)
			}
			if pos, err := strconv.Atoi(field[1]); err == nil {
				for i, alt := range alts {
					var snv Snv
					if alt == "*" {
						continue
					}
					snv.Variant = core.Variant{
						Chrom: strings.Replace(field[0], "chr", "", 1),
						Start: pos,
						End:   0,
						Ref:   core.Sequence(field[3]),
						Alt:   core.Sequence(alt),
					}
					snv.Variant.ConvertSnv()
					snv.Information = Information{
						Depth:      depth,
						Qual:       qual,
						GatkFilter: gatkFilter,
						Genotype:   genotypes[i],
						Ratio:      ratios[i],
					}
					snv.OtherInfo = string(line)
					snvs = append(snvs, snv)
				}
			}
		}
	}
	sort.Sort(snvs)
	return snvs
}
