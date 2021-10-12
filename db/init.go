package db

import (
	"fmt"
	"grandanno/seq"
	"log"
)

var DB Database

var ChromArray Chromosomes

var UpDownStreamLen int

var RefIndexStepLen int

var SymbolToEntrezId EntrezIdMap

func InitDatabase(databasePath string, genomeVersion string) {
	DB = NewDatabase(databasePath, genomeVersion)
	DB.Check()
	ChromArray = NewChromosomesFromReferenceDict(DB.ReferenceDictFile())
}

func PrepareDatabase(databasePath string, genomeVersion string) {
	DB = NewDatabase(databasePath, genomeVersion)
	ChromArray = NewChromosomesFromReferenceDict(DB.ReferenceDictFile())
	// mRNA
	log.Printf("start read %s\n", DB.ReferenceFastaFile())
	reference := seq.NewFasta(DB.ReferenceFastaFile())
	fmt.Printf("start read %s\n", DB.RefgeneFile())
	refgenes := NewRefgenes(DB.RefgeneFile())
	fmt.Printf("start write %s\n", DB.MrnaFastaFile())
	WriteMranFasta(refgenes, reference, DB.MrnaFastaFile())
	// Reference Index
	fmt.Printf("start init and write %s\n", DB.RefIndexFile())
	CreateRefIndex(refgenes, DB.RefIndexFile())
}

func init() {
	UpDownStreamLen = 3000
	RefIndexStepLen = 300000
}
