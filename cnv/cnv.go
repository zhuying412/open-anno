package cnv

import "grandanno/core"

type Cnv interface {
	GetVariant() core.Variant
	GetTypo() string
}

type Cnvs []Cnv

func (cnvs Cnvs) Len() int {
	return len(cnvs)
}

func (cnvs Cnvs) Less(i, j int) bool {
	starti, endi := cnvs[i].GetVariant().GetDigitalPosition()
	startj, endj := cnvs[j].GetVariant().GetDigitalPosition()
	if starti == startj {
		return endi < endj
	} else {
		return starti < startj
	}
}

func (cnvs Cnvs) Swap(i, j int) {
	cnvs[i], cnvs[j] = cnvs[j], cnvs[i]
}
