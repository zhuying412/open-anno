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

type GatkSnv struct {
	Variant     core.Variant
	Information Information
	OtherInfo   string
}

func (snv GatkSnv) GetVariant() core.Variant {
	return snv.Variant
}

func (snv GatkSnv) GetTypo() string {
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

type GatkVcfLine struct {
	Chrom        string
	Pos          int
	Ref          string
	Alts         []string
	GatkFilter   string
	Qual         float32
	InfoList     []string
	FormatKeys   []string
	FormatValues []string
	Line         string
}

func (vcfLine *GatkVcfLine) readLine(line string) error {
	field := strings.Split(line, "\t")
	vcfLine.Line = line
	vcfLine.Chrom = strings.Replace(field[0], "chr", "", 1)
	vcfLine.Alts = strings.Split(field[4], ",")
	vcfLine.Ref = field[3]
	vcfLine.GatkFilter = field[6]
	vcfLine.InfoList = strings.Split(field[7], ";")
	vcfLine.FormatKeys = strings.Split(field[len(field)-2], ":")
	vcfLine.FormatValues = strings.Split(field[len(field)-1], ":")
	if qual, err := strconv.ParseFloat(field[5], 32); err == nil {
		vcfLine.Qual = float32(qual)
	} else {
		vcfLine.Qual = float32(-1)
	}
	if pos, err := strconv.Atoi(field[1]); err == nil {
		vcfLine.Pos = pos
	} else {
		return err
	}
	return nil
}

func (vcfLine GatkVcfLine) getGenotypes() []float32 {
	count := len(vcfLine.Alts)
	genotypes := make([]float32, count)
	for i := 0; i < count; i++ {
		genotypes[i] = float32(-1)
	}
	for _, info := range vcfLine.InfoList {
		if strings.HasPrefix(info, "AF=") {
			for i, gt := range strings.Split(info[3:], ",") {
				if i > count {
					break
				}
				if _gt, err := strconv.ParseFloat(gt, 32); err == nil {
					genotypes[i] = float32(_gt)
				}
			}
			break
		}
	}
	return genotypes
}

func (vcfLine GatkVcfLine) getDepthRatios() (int, []float32) {
	count := len(vcfLine.Alts)
	ratios := make([]float32, count)
	depth := -1
	for i := 0; i < count; i++ {
		ratios[i] = float32(-1)
	}
	for i, key := range vcfLine.FormatKeys {
		if key == "DP" {
			if dp, err := strconv.Atoi(vcfLine.FormatValues[i]); err == nil {
				depth = dp
			}
		}
		if key == "AD" {
			sum := 0
			counts := make([]int, count)
			_counts := strings.Split(vcfLine.FormatValues[i], ",")
			for i := 0; i <= count && i < len(_counts); i++ {
				if c, err := strconv.Atoi(_counts[i]); err == nil {
					sum += c
					if i > 0 {
						counts[i-1] = c
					}
				}
			}
			if sum > 0 {
				for i := 0; i < count; i++ {
					ratios[i] = float32(counts[i]) / float32(sum)
				}
			}
		}
	}
	return depth, ratios
}

func (vcfLine GatkVcfLine) getSnvs() Snvs {
	var snvs Snvs
	genotypes := vcfLine.getGenotypes()
	depth, ratios := vcfLine.getDepthRatios()
	for i, alt := range vcfLine.Alts {
		if alt == "*" {
			continue
		}
		snv := GatkSnv{
			Variant: core.Variant{
				Chrom: vcfLine.Chrom,
				Start: vcfLine.Pos,
				End:   0,
				Ref:   core.Sequence(vcfLine.Ref),
				Alt:   core.Sequence(alt),
			},
			Information: Information{
				Depth:      depth,
				Qual:       vcfLine.Qual,
				GatkFilter: vcfLine.GatkFilter,
				Genotype:   genotypes[i],
				Ratio:      ratios[i],
			},
			OtherInfo: vcfLine.Line,
		}
		snv.Variant.ConvertSnv()
		snvs = append(snvs, snv)
	}
	return snvs
}

func (snvs *Snvs) ReadGatkVcfFile(vcfFile string) {
	var gatkVcfLine GatkVcfLine
	if fp, err := os.Open(vcfFile); err == nil {
		defer fp.Close()
		var lines [][]byte
		if strings.HasSuffix(strings.ToLower(vcfFile), ".gz") {
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
			if err := gatkVcfLine.readLine(string(line)); err == nil {
				*snvs = append(*snvs, gatkVcfLine.getSnvs()...)
			}
		}
	} else {
		panic(err)
	}
	sort.Sort(*snvs)
}
