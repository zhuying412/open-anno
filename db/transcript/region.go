package transcript

import "strings"

type Region struct {
	Start     int    `json:"start"`
	End       int    `json:"end"`
	Type      string `json:"type"`
	ExonOrder int    `json:"exon_order"`
}

func (r Region) IsCDS() bool {
	return r.Type == "CDS"
}

func (r Region) IsUTR() bool {
	return strings.HasPrefix(r.Type, "UTR")
}

func (r Region) IsIntron() bool {
	return r.Type == "intron"
}

type Regions []Region

func (r Regions) GetPrev(currentIndex int, strand byte) (Region, bool) {
	var index int
	if strand == '+' {
		index = currentIndex - 1
	} else {
		index = currentIndex + 1
	}
	if index >= 0 && index < r.Len() {
		return r[index], true
	}
	return Region{}, false
}

func (r Regions) GetNext(currentIndex int, strand byte) (Region, bool) {
	var index int
	if strand == '+' {
		index = currentIndex + 1
	} else {
		index = currentIndex - 1
	}
	if index >= 0 && index < r.Len() {
		return r[index], true
	}
	return Region{}, false
}

func (r Regions) FindOne(pos int, strand byte) (index int, exonLen int) {
	for j := 0; j < r.Len(); j++ {
		index = j
		if strand == '-' {
			index = r.Len() - j - 1
		}
		region := r[index]
		if region.Start <= pos && pos <= region.End {
			break
		}
		if region.IsCDS() {
			exonLen += region.End - region.Start + 1
		}
	}
	return index, exonLen
}

func (r Regions) FindMany(start int, end int, strand byte) (indexes []int, lenL int, lenR int) {
	for j := 0; j < r.Len(); j++ {
		i := j
		if strand == '-' {
			i = r.Len() - j - 1
		}
		region := r[i]
		condition1, condition2 := region.Start > end, region.End < start
		if strand == '-' {
			condition1, condition2 = region.End < start, region.Start > end
		}
		if condition1 {
			if region.IsCDS() {
				lenR += region.End - region.Start + 1
			}
		} else if condition2 {
			if region.IsCDS() {
				lenL += region.End - region.Start + 1
			}
		} else {
			indexes = append(indexes, i)
		}
	}
	return indexes, lenL, lenR
}

func (r Regions) HasCds() bool {
	for _, region := range r {
		if region.IsCDS() {
			return true
		}
	}
	return false
}

func (r Regions) Len() int {
	return len(r)
}

func (r Regions) Less(i, j int) bool {
	return r[i].Start < r[j].Start
}

func (r Regions) Swap(i, j int) {
	r[i], r[j] = r[j], r[i]
}
