package scheme

import (
	"fmt"
	"sort"
)

const (
	RType_CDS    = "CDS"
	RType_INTRON = "intron"
	RType_UTR    = "UTR"
)

var IS_EXON_REGION = false

type Region struct {
	Chrom    string `json:"chrom"`
	Start    int    `json:"start"`
	End      int    `json:"end"`
	Type     string `json:"type"`
	Order    int    `json:"order"`
	Exon     string `json:"exon"`
	Sequence string `json:"sequence"`
}

func (this *Region) Name() string {
	if !IS_EXON_REGION {
		if this.Order == 0 {
			return this.Type
		}
		return fmt.Sprintf("%s%d", this.Type, this.Order)
	}
	return this.Exon

}

func (this Region) Equal(r Region) bool {
	return this.Chrom == r.Chrom && this.Start == r.Start && this.End == r.End
}

func (this Region) Exists() bool {
	return this.Chrom != "" && this.Start != 0 && this.End != 0
}

type Regions []Region

func (this Regions) Len() int           { return len(this) }
func (this Regions) Swap(i, j int)      { this[i], this[j] = this[j], this[i] }
func (this Regions) Less(i, j int) bool { return this[i].Start < this[j].Start }

func NewRegions(trans Transcript) (Regions, error) {
	regions := make(Regions, 0)
	for i := 0; i < trans.ExonCount; i++ {
		exon := fmt.Sprintf("exon%d", i+1)
		if trans.Strand == "-" {
			exon = fmt.Sprintf("exon%d", trans.ExonCount-i)
		}
		if i > 0 && trans.ExonEnds[i-1]+1 <= trans.ExonStarts[i]-1 {
			regions = append(regions, Region{
				Chrom: trans.Chrom,
				Start: trans.ExonEnds[i-1] + 1,
				End:   trans.ExonStarts[i] - 1,
				Type:  RType_INTRON,
			})
		}
		start, end := trans.ExonStarts[i], trans.ExonEnds[i]
		if trans.CdsStart > end || trans.CdsEnd < start {
			regions = append(regions, Region{
				Chrom: trans.Chrom,
				Start: start,
				End:   end,
				Type:  RType_UTR,
				Exon:  exon,
			})
		} else if trans.CdsStart <= start && end <= trans.CdsEnd {
			regions = append(regions, Region{
				Chrom: trans.Chrom,
				Start: start,
				End:   end,
				Type:  RType_CDS,
				Exon:  exon,
			})
		} else {
			newStart, newEnd := start, end
			hasUTR1, hasUTR2 := false, false
			if start < trans.CdsStart && trans.CdsStart < end {
				newStart = trans.CdsStart
				hasUTR1 = true
			}
			if start < trans.CdsEnd && trans.CdsEnd < end {
				newEnd = trans.CdsEnd
				hasUTR2 = true
			}
			if hasUTR1 {
				regions = append(regions, Region{
					Chrom: trans.Chrom,
					Start: start,
					End:   trans.CdsStart - 1,
					Type:  RType_UTR,
					Exon:  exon,
				})
			}
			regions = append(regions, Region{
				Chrom: trans.Chrom,
				Start: newStart,
				End:   newEnd,
				Type:  RType_CDS,
				Exon:  exon,
			})
			if hasUTR2 {
				regions = append(regions, Region{
					Chrom: trans.Chrom,
					Start: trans.CdsEnd + 1,
					End:   end,
					Type:  RType_UTR,
					Exon:  exon,
				})
			}
		}
	}
	if trans.Strand == "+" {
		sort.Sort(regions)
	} else {
		sort.Sort(sort.Reverse(regions))
	}
	intron, cds, utr := 1, 1, 5
	var err error
	for i := 0; i < len(regions); i++ {
		switch regions[i].Type {
		case RType_INTRON:
			regions[i].Order = intron
			intron++
		case RType_CDS:
			regions[i].Order = cds
			cds++
			utr = 3
		case RType_UTR:
			regions[i].Order = utr
		}
	}
	sort.Sort(regions)
	return regions, err
}
