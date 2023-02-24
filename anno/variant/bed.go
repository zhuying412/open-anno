package variant

import (
	"open-anno/pkg"
	"strconv"
	"strings"
)

type CNVs []AnnoVariant

func (this CNVs) Len() int           { return len(this) }
func (this CNVs) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this CNVs) Less(i, j int) bool { return this[i].Less(this[j]) }

func (this CNVs) AggregateByChrom() (map[string]CNVs, error) {
	cnvs := make(map[string]CNVs)
	for _, cnv := range this {
		chrom := cnv.Chrom
		if rows, ok := cnvs[chrom]; ok {
			cnvs[chrom] = append(rows, cnv)
		} else {
			cnvs[chrom] = CNVs{cnv}
		}
	}
	return cnvs, nil
}

func ReadBED(bedFile string) (map[string]CNVs, error) {
	cnvs := make(map[string]CNVs)
	reader, err := pkg.NewIOReader(bedFile)
	if err != nil {
		return cnvs, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		start, err := strconv.Atoi(row[1])
		if err != nil {
			return cnvs, err
		}
		end, err := strconv.Atoi(row[2])
		if err != nil {
			return cnvs, err
		}
		cnv := AnnoVariant{
			Chrom: row[0],
			Start: start,
			End:   end,
			Ref:   "DIP",
			Alt:   row[3],
		}
		if rows, ok := cnvs[cnv.Chrom]; ok {
			cnvs[cnv.Chrom] = append(rows, cnv)
		} else {
			cnvs[cnv.Chrom] = CNVs{cnv}
		}
	}
	return cnvs, nil
}
