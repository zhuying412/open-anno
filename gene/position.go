package gene

import (
	"strconv"
	"strings"
)

type Position struct {
	ExonStart  int   `json:"exon_start"`
	ExonEnd    int   `json:"exon_end"`
	CdsStart   int   `json:"cds_start"`
	CdsEnd     int   `json:"cds_end"`
	ExonStarts []int `json:"exon_starts"`
	ExonEnds   []int `json:"exon_ends"`
}

func NewPosition(refgeneLineField []string) Position {
	position := Position{}
	var err error
	if position.ExonStart, err = strconv.Atoi(refgeneLineField[4]); err != nil {
		panic(err)
	} else {
		position.ExonStart += 1
	}
	if position.ExonEnd, err = strconv.Atoi(refgeneLineField[5]); err != nil {
		panic(err)
	}
	if position.CdsStart, err = strconv.Atoi(refgeneLineField[6]); err != nil {
		panic(err)
	} else {
		position.CdsStart += 1
	}
	if position.CdsEnd, err = strconv.Atoi(refgeneLineField[7]); err != nil {
		panic(err)
	}
	for _, strPos := range strings.Split(strings.Trim(refgeneLineField[9], ","), ",") {
		if pos, err := strconv.Atoi(strPos); err == nil {
			position.ExonStarts = append(position.ExonStarts, pos+1)
		}
	}
	for _, strPos := range strings.Split(strings.Trim(refgeneLineField[10], ","), ",") {
		if pos, err := strconv.Atoi(strPos); err == nil {
			position.ExonEnds = append(position.ExonEnds, pos)
		}
	}
	return position
}
