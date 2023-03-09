package db

import (
	"io"
	"open-anno/anno"
	"open-anno/pkg"
	"strconv"
	"strings"

	"github.com/brentp/vcfgo"
)

type RegionVar struct {
	Chrom string `json:"chrom"`
	Start int    `json:"start"`
	End   int    `json:"end"`
	Info  string `json:"name"`
}

type RegionVarScanner struct {
	pkg.IOScanner
	Name string
}

func NewRegionVarScanner(reader io.ReadCloser) RegionVarScanner {
	scanner := pkg.NewIOScanner(reader)
	scanner.Scan()
	return RegionVarScanner{IOScanner: scanner, Name: strings.Split(scanner.Text(), "\t")[3]}
}

func (this RegionVarScanner) VCFHeaderInfo() map[string]*vcfgo.Info {
	return map[string]*vcfgo.Info{
		this.Name: {
			Id:          this.Name,
			Description: "Annotation of " + this.Name,
			Number:      ".",
			Type:        "UNKNOWN",
		},
	}
}

func (this RegionVarScanner) Row() (RegionVar, error) {
	fields := strings.Split(this.Text(), "\t")
	start, err := strconv.Atoi(fields[1])
	if err != nil {
		return RegionVar{}, err
	}
	end, err := strconv.Atoi(fields[2])

	return RegionVar{Chrom: fields[0], Start: start, End: end, Info: fields[3]}, err
}

// ReadAll 读取RegionBased Bed
func (this RegionVarScanner) ReadAll() (map[string][]RegionVar, error) {
	regVars := make(map[string][]RegionVar, 0)
	for this.Scan() {
		row, err := this.Row()
		if err != nil {
			return regVars, err
		}
		if vars, ok := regVars[row.Chrom]; ok {
			regVars[row.Chrom] = append(vars, row)
		} else {
			regVars[row.Chrom] = []RegionVar{row}
		}
	}
	return regVars, nil
}

// AnnoRegionBased 注释RegionBased
func AnnoRegionBased(variants anno.Variants, regVars map[string][]RegionVar, name string, overlap float64) (anno.AnnoInfos, error) {
	annoInfos := make(anno.AnnoInfos)
	// 读取DBFile DB文件
	// log.Printf("Read DBFile: %s ...", dbFile)
	// 开始注释
	for _, v := range variants {
		annoVariant := v.AnnoVariant()
		if regs, ok := regVars[annoVariant.Chrom]; ok {
			var infos []string
			for _, reg := range regs {
				if annoVariant.End >= reg.Start && annoVariant.Start <= reg.End {
					vlen := annoVariant.End - annoVariant.Start + 1
					olen := pkg.Min(annoVariant.End, reg.End) - pkg.Max(annoVariant.Start, reg.Start) + 1
					if float64(olen)/float64(vlen) >= overlap {
						infos = append(infos, reg.Info)
					}
				}
			}
			if len(infos) > 0 {
				annoInfos[annoVariant.PK()] = []anno.AnnoInfo{{Key: name, Value: strings.Join(infos, ",")}}
			}
		}
	}
	return annoInfos, nil
}
