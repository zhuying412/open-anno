package db

import (
	"open-anno/pkg"
	"strconv"
	"strings"
)

// RegionBed RegionBased Bed
type RegionVar struct {
	Chrom string `json:"chrom"`
	Start int    `json:"start"`
	End   int    `json:"end"`
	Info  string `json:"info"`
}

// ReadRegionBed 读取RegionBased Bed
func ReadRegionBased(dbFile string) (map[string][]RegionVar, string, error) {
	var regMap map[string][]RegionVar
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
		reg := RegionVar{Chrom: row[0], Info: row[3]}
		reg.Start, err = strconv.Atoi(row[1])
		if err != nil {
			return regMap, name, err
		}
		reg.End, err = strconv.Atoi(row[2])
		if err != nil {
			return regMap, name, err
		}
		if _, ok := regMap[reg.Chrom]; ok {
			regMap[reg.Chrom] = make([]RegionVar, 0)
		}
		regMap[reg.Chrom] = append(regMap[reg.Chrom], reg)
	}
	return regMap, name, nil
}
