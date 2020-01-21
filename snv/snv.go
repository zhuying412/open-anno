package snv

import "grandanno/core"

type Snv interface {
	GetVariant() core.Variant
	GetTypo() string
}

type Snvs []Snv

func (snvs Snvs) Len() int {
	return len(snvs)
}

func (snvs Snvs) Less(i, j int) bool {
	starti, endi := snvs[i].GetVariant().GetDigitalPosition()
	startj, endj := snvs[j].GetVariant().GetDigitalPosition()
	if starti == startj {
		return endi < endj
	} else {
		return starti < startj
	}
}

func (snvs Snvs) Swap(i, j int) {
	snvs[i], snvs[j] = snvs[j], snvs[i]
}
