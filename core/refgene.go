package core

import (
	"fmt"
)

type Position struct {
	ExonStart  int
	ExonEnd    int
	CdsStart   int
	CdsEnd     int
	ExonStarts []int
	ExonEnds   []int
}

type Region struct {
	Start     int
	End       int
	Typo      string
	ExonOrder int
}

type Regions []Region

type Refgene struct {
	Chrom      string
	Strand     bool
	Gene       string
	EntrezId   string
	Transcript string
	Position   Position
	Regions    Regions
	Tag        string
}

type Refgenes []Refgene

func (refgene Refgene) GetSn() string {
	return fmt.Sprintf("%s|%s:%d:%d", refgene.Transcript, refgene.Chrom, refgene.Position.ExonStart, refgene.Position.ExonEnd)
}

func (refgenes *Refgenes) Read() {

}
