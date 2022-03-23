package regionbased

import (
	"bufio"
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/variant"
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

func ReadDB(dbfile string) (RegionBaseds, string, error) {
	regionBaseds := make(RegionBaseds, 0)
	fi, err := os.Open(dbfile)
	if err != nil {
		log.Fatal(err)
	}
	defer fi.Close()
	scanner := bufio.NewScanner(fi)
	var header string
	for scanner.Scan() {
		line := scanner.Text()
		if header == "" {
			header = line
			continue
		}
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
			Chrom:     fields[0],
			Start:     start,
			End:       end,
			Otherinfo: fields[3],
		})
	}
	sort.Sort(regionBaseds)
	return regionBaseds, header, err
}

func Anno(variants variant.Variants, dbfile string, overlap float64, writer *os.File) {
	sort.Sort(variants)
	regionBaseds, header, err := ReadDB(dbfile)
	if err != nil {
		log.Fatal(err)
	}
	writer.WriteString(header + "\n")
	for _, variant := range variants {
		var info []string
		for _, dbvar := range regionBaseds {
			if variant.End >= dbvar.Start && variant.Start <= dbvar.End {
				vlen := variant.End - variant.Start + 1
				olen := pkg.Max(variant.Start, dbvar.Start) - pkg.Min(variant.Start, dbvar.Start) + 1
				if float64(olen)/float64(vlen) >= overlap {
					info = append(info, dbvar.Otherinfo)
				}
			}
		}
		anno := strings.Join(info, ",")
		if anno == "" {
			anno = "."
		}
		writer.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s\t%s\t%s\n",
			variant.Chrom, variant.Start, variant.End,
			variant.Ref, variant.Alt, anno,
		))
	}
}
