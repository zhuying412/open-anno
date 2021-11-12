package chromosome

import (
	"log"
)

type Chromosome struct {
	Name   string
	Length int
}

type Chromosomes []Chromosome

func (c Chromosomes) GetByName(name string) (int, Chromosome) {
	for order, chrom := range c {
		if chrom.Name == name {
			return order, chrom
		}
	}
	log.Panicf("Not found chromosome:%s", name)
	return 0, Chromosome{}
}
