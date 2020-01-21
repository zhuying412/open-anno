package prepare

import (
	"fmt"
	"grandanno/core"
	"os"
	"sort"
	"strings"
)

type Refidx struct {
	Chrom       string
	Start       int
	End         int
	Transcripts []string
}

type Refidxs []Refidx

type RefidxDict map[string]Refidxs

func (refidx Refidx) GetDigitalPosition() (int, int) {
	start := core.ChromOrderDict[refidx.Chrom]*1e9 + refidx.Start
	end := core.ChromOrderDict[refidx.Chrom]*1e9 + refidx.End
	return start, end
}

func (refidxs Refidxs) Len() int {
	return len(refidxs)
}

func (refidxs Refidxs) Less(i, j int) bool {
	starti, endi := refidxs[i].GetDigitalPosition()
	startj, endj := refidxs[j].GetDigitalPosition()
	if starti == startj {
		return endi < endj
	} else {
		return starti < startj
	}
}

func (refidxs Refidxs) Swap(i, j int) {
	refidxs[i], refidxs[j] = refidxs[j], refidxs[i]
}

func (refidxs Refidxs) SetTranscript(refgenes Refgenes, refidxsChannel chan Refidxs) {
	for i := 0; i < len(refidxs); i++ {
		for j := 0; j < len(refgenes); j++ {
			starti, endi := refidxs[i].GetDigitalPosition()
			startj, endj := refgenes[j].GetDigitalPosition()
			if starti <= endj && endi >= startj {
				refidxs[i].Transcripts = append(refidxs[i].Transcripts, refgenes[j].GetSn())
			}
		}
	}
	refidxsChannel <- refidxs
}

func (refidxDict RefidxDict) Init(refidxStepLen int) {
	for chrom, length := range core.ChromLengthDict {
		for i := 0; i < length; i += refidxStepLen {
			end := i + refidxStepLen
			if end > length {
				end = length
			}
			refdix := Refidx{Chrom: chrom, Start: i + 1, End: end}
			if refidxs, ok := refidxDict[chrom]; ok {
				refidxDict[chrom] = append(refidxs, refdix)
			} else {
				refidxDict[chrom] = Refidxs{refdix}
			}
		}
	}
}

func (refidxDict RefidxDict) Write(refidxFile string, refgeneDict RefgeneDict) {
	fo, err := os.Create(refidxFile)
	if err == nil {
		defer fo.Close()
		refidxsChannel := make(chan Refidxs)
		for _, chrom := range core.ChromList {
			refidxs := refidxDict[chrom]
			refgenes := refgeneDict[chrom]
			sort.Sort(refidxs)
			sort.Sort(refgenes)
			go refidxs.SetTranscript(refgenes, refidxsChannel)
		}
		for range core.ChromList {
			if refidxs, ok := <-refidxsChannel; ok {
				refidxDict[refidxs[0].Chrom] = refidxs
			} else {
				break
			}
		}
		for _, chrom := range core.ChromList {
			for _, refidx := range refidxDict[chrom] {
				if len(refidx.Transcripts) > 0 {
					if _, err := fo.WriteString(fmt.Sprintf(
						"%s\t%d\t%d\t%s\n",
						refidx.Chrom, refidx.Start, refidx.End,
						strings.Join(refidx.Transcripts, ","),
					)); err != nil {
						panic(err.Error())
					}
				}
			}
		}
		close(refidxsChannel)
	}
}
