package snv

import (
	"bufio"
	"bytes"
	"grandanno/core"
	"io"
	"os"
	"strconv"
	"strings"
)

type Vcf struct {
	Info   map[string][]string
	Format map[string][]string
}

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

func (vcf Vcf) readInfo(info string) {
	for _, s := range strings.Split(info, ";") {
		if strings.Contains(s, "=") {
			arr := strings.Split(s, "=")
			vcf.Info[arr[0]] = strings.Split(arr[1], ",")
		}
	}
}

func (vcf Vcf) readFormat(formatKey string, formatValue string) {
	keys := strings.Split(formatKey, ":")
	values := strings.Split(formatValue, ":")
	for i, key := range keys {
		vcf.Format[key] = strings.Split(values[i], ",")
	}
}

func (vcf Vcf) getGenotypes(n int) []float32 {
	genotypes := make([]float32, n)
	if values, ok := vcf.Info["AF"]; ok {
		for i, value := range values {
			if valueFloat, err := strconv.ParseFloat(value, 32); err == nil {
				genotypes[i] = float32(valueFloat)
			}
		}
	} else {
		for i := 0; i < n; i++ {
			genotypes[i] = -1
		}
	}
	return genotypes
}

func (vcf Vcf) getDepthRatio(n int) (int, []float32) {
	ratios := make([]float32, n)
	depth := -1
	depthList, ok1 := vcf.Format["DP"]
	countList, ok2 := vcf.Format["AD"]
	if ok1 && ok2 && len(depthList) > 0 && len(countList) > n {
		if d, err := strconv.Atoi(depthList[0]); err == nil {
			depth = d
		}
		sum := 0
		countIntList := make([]int, n)
		for i := 0; i <= n; i += n {
			if c, err := strconv.Atoi(countList[i]); err == nil {
				sum += c
				if i > 0 {
					countIntList[i-1] = c
				}
			}
		}
		if sum > 0 {
			for i := 1; i < n; i += n {
				ratios[i] = float32(countIntList[i] / sum)
			}
		}
	} else {
		for i := 0; i < n; i++ {
			ratios[i] = -1
		}
	}
	return depth, ratios
}

func (vcf Vcf) ReadLine(vcfLine string) Snvs {
	field := strings.Split(vcfLine, "\t")
	var alts []string
	for _, alt := range strings.Split(field[4], ",") {
		alts = append(alts, alt)
	}
	vcf.readInfo(field[7])
	vcf.readFormat(field[len(field)-1], field[len(field)])
	genotypes := vcf.getGenotypes(len(alts))
	depth, ratios := vcf.getDepthRatio(len(alts))
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

func (vcf Vcf) Read(vcfFile string) Snvs {
	var snvs Snvs
	if fp, err := os.Open(vcfFile); err == nil {
		defer fp.Close()
		reader := bufio.NewReader(fp)
		for {
			if line, err := reader.ReadBytes('\n'); err == nil {
				line = bytes.TrimSpace(line)
				if len(line) == 0 {
					continue
				}
				tmpSnvs := vcf.ReadLine(string(line))
				snvs = append(snvs, tmpSnvs...)
			} else {
				if err == io.EOF {
					break
				} else {
					panic(err.Error())
				}
			}
		}
	} else {
		panic(err.Error())
	}
	return snvs
}

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
