package variant

import (
	"bufio"
	"log"
	"open-anno/pkg"
	"os"
	"sort"
	"strconv"
	"strings"
)

type RegionBased struct {
	Chrom     string `json:"chrom"`
	Start     int    `json:"start"`
	End       int    `json:"end"`
	Otherinfo string `json:"otherinfo"`
}

type RegionBaseds []RegionBased

func (this RegionBaseds) Len() int      { return len(this) }
func (this RegionBaseds) Swap(i, j int) { this[i], this[j] = this[j], this[i] }
func (this RegionBaseds) Less(i, j int) bool {
	if this[i].Start == this[j].Start {
		return this[i].End < this[j].End
	} else {
		return this[i].Start < this[j].Start
	}
}

func ReadRegionBasedLine(line string) (RegionBased, error) {
	fields := strings.Split(line, "\t")
	variant := RegionBased{
		Chrom:     pkg.FormatChrom(fields[0]),
		Otherinfo: fields[3],
	}
	var err error
	variant.Start, err = strconv.Atoi(fields[1])
	if err != nil {
		return variant, err
	}
	variant.End, err = strconv.Atoi(fields[2])
	if err != nil {
		return variant, err
	}
	return variant, err
}

func ReadRegionBasedDB(infile string) (RegionBaseds, string, error) {
	regionBaseds := make(RegionBaseds, 0)
	fi, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer fi.Close()
	scanner := bufio.NewScanner(fi)
	scanner.Scan()
	header := scanner.Text()
	if !strings.HasPrefix(header, "#Chr") {
		log.Fatalf("error database file, header not found: %s", infile)
	}
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return regionBaseds, header, err
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return regionBaseds, header, err
		}
		regionBaseds = append(regionBaseds, RegionBased{
			Chrom:     pkg.FormatChrom(fields[0]),
			Start:     start,
			End:       end,
			Otherinfo: fields[3],
		})
	}
	sort.Sort(regionBaseds)
	return regionBaseds, header, err
}
