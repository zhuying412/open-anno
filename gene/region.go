package gene

type Region struct {
	Start     int    `json:"start"`
	End       int    `json:"end"`
	Type      string `json:"type"`
	ExonOrder int    `json:"exon_order"`
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

func (r Regions) Len() int {
	return len(r)
}

func (r Regions) Less(i, j int) bool {
	return r[i].Start < r[j].Start
}

func (r Regions) Swap(i, j int) {
	r[i], r[j] = r[j], r[i]
}
