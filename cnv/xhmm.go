package cnv

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"grandanno/core"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

type Information struct {
	MeanDepthDepth    float32
	MeanOriginalDepth float32
}

type XhmmCnv struct {
	Variant     core.Variant
	Information Information
	OtherInfo   string
}

func (cnv XhmmCnv) GetVariant() core.Variant {
	return cnv.Variant
}

func (cnv XhmmCnv) GetTypo() string {
	return strings.Trim(string(cnv.Variant.Alt), "<>")
}

type XhmmCnvDict map[string]Cnvs

type XhmmVcfLine struct {
	Chrom     string
	Start     int
	End       int
	Ref       string
	Alts      []string
	OtherInfo []string
}

func (info Information) read(otherInfo string) int {
	tmp := strings.Split(otherInfo, ":")
	lenght := len(tmp)
	if mrd, err := strconv.ParseFloat(tmp[lenght-3], 32); err == nil {
		info.MeanDepthDepth = float32(mrd)
	} else {
		info.MeanDepthDepth = float32(-1)
	}
	if mord, err := strconv.ParseFloat(tmp[lenght-1], 32); err == nil {
		info.MeanOriginalDepth = float32(mord)
	} else {
		info.MeanDepthDepth = float32(-1)
	}
	if genotype, err := strconv.Atoi(tmp[0]); err == nil {
		return genotype
	}
	return 0
}

func (xhmmVcfLine *XhmmVcfLine) getCnvDict(head []string) map[string]Cnv {
	cnvDict := make(map[string]Cnv)
	for i, otherInfo := range xhmmVcfLine.OtherInfo {
		var info Information
		sample := head[i]
		if genotype := info.read(otherInfo); genotype != 0 {
			cnv := XhmmCnv{
				Variant: core.Variant{
					Chrom: xhmmVcfLine.Chrom,
					Start: xhmmVcfLine.Start,
					End:   xhmmVcfLine.End,
					Ref:   core.Sequence(xhmmVcfLine.Ref),
					Alt:   core.Sequence(xhmmVcfLine.Alts[genotype-1]),
				},
				Information: info,
				OtherInfo:   "",
			}
			cnvDict[sample] = cnv
		}
	}
	return cnvDict
}

func (xhmmVcfLine *XhmmVcfLine) readLine(vcfLine string) error {
	field := strings.Split(vcfLine, "\t")
	tmp := strings.Split(field[2], ":")
	xhmmVcfLine.Chrom = tmp[0]
	tmp = strings.Split(tmp[1], "-")
	if start, err := strconv.Atoi(tmp[0]); err == nil {
		xhmmVcfLine.Start = start
		if end, err := strconv.Atoi(tmp[1]); err == nil {
			xhmmVcfLine.End = end
			xhmmVcfLine.Ref = field[3]
			xhmmVcfLine.Alts = strings.Split(field[4], ",")
			xhmmVcfLine.OtherInfo = field[9:]
		} else {
			return err
		}
	} else {
		return err
	}
	return nil
}

func (xhmmCnvDict XhmmCnvDict) ReadXhmmVcfFile(vcfFile string) {
	var head []string
	var xhmmVcfLine XhmmVcfLine
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
			if len(line) == 0 {
				continue
			}
			if line[0] == '#' {
				if bytes.HasPrefix(line, []byte("#CHROM")) {
					head = strings.Split(string(line), "\t")[9:]
				}
			} else {
				if err := xhmmVcfLine.readLine(string(line)); err != nil {
					cnvDict := xhmmVcfLine.getCnvDict(head)
					for sample, cnv := range cnvDict {
						if cnvs, ok := xhmmCnvDict[sample]; ok {
							xhmmCnvDict[sample] = append(cnvs, cnv)
						} else {
							xhmmCnvDict[sample] = Cnvs{cnv}
						}
					}
				} else {
					panic(err)
				}

			}
		}
	}
}
