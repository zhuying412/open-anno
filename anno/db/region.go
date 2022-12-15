package db

import (
	"fmt"
	"log"
	"open-anno/anno"
	"open-anno/pkg"
	"strconv"
	"strings"
)

// RegionBed RegionBased Bed
type RegionBed struct {
	Chrom string `json:"chrom"`
	Start int    `json:"start"`
	End   int    `json:"end"`
	Info  string `json:"Info"`
}

// ReadRegionBed 读取RegionBased Bed
func ReadRegionBed(dbFile string) (map[string][]RegionBed, string, error) {
	var regMap map[string][]RegionBed
	var name string
	reader, err := pkg.NewIOReader(dbFile)
	if err != nil {
		return regMap, name, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	scanner.Scan()
	name = strings.Split(scanner.Text(), "\t")[3]
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		reg := RegionBed{Chrom: row[0], Info: row[3]}
		reg.Start, err = strconv.Atoi(row[1])
		if err != nil {
			return regMap, name, err
		}
		reg.End, err = strconv.Atoi(row[2])
		if err != nil {
			return regMap, name, err
		}
		if _, ok := regMap[reg.Chrom]; ok {
			regMap[reg.Chrom] = make([]RegionBed, 0)
		}
		regMap[reg.Chrom] = append(regMap[reg.Chrom], reg)
	}
	return regMap, name, nil
}

// AnnoRegionBased 注释RegionBased
func AnnoRegionBased(annoInputFile, dbFile, annoOutputFile string, overlap float64) error {
	// 读取变异输入文件
	log.Printf("Read AnnoInput: %s ...", annoInputFile)
	variantMap, err := anno.ReadAnnoInput(annoInputFile)
	if err != nil {
		return err
	}
	// 读取DBFile DB文件
	log.Printf("Read DBFile: %s ...", dbFile)
	regMap, name, err := ReadRegionBed(dbFile)
	if err != nil {
		return err
	}
	// 定义输出
	writer, err := pkg.NewIOWriter(annoInputFile)
	if err != nil {
		return err
	}
	defer writer.Close()
	// 开始注释
	// 输出表头
	fmt.Fprintf(writer, "Chr\tStart\tEnd\tRef\tAlt\t%s\n", name)
	for chrom, variants := range variantMap {
		if regs, ok := regMap[chrom]; ok {
			for _, variant := range variants {
				var annos []string
				for _, reg := range regs {
					if variant.End >= reg.Start && variant.Start <= reg.End {
						vlen := variant.End - variant.Start + 1
						olen := pkg.Min(variant.End, reg.End) - pkg.Max(variant.Start, reg.Start) + 1
						if float64(olen)/float64(vlen) >= overlap {
							annos = append(annos, reg.Info)
						}
					}
				}
				if len(annos) > 0 {
					fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\n",
						variant.Chrom, variant.Start, variant.End,
						variant.Ref, variant.Alt, strings.Join(annos, ","))
				}
			}
		}
	}
	return err
}
