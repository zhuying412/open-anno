package db

import (
	"github.com/spf13/viper"
	"grandanno/seq"
	"log"
	"os"
	"path"
)

var ChromArray Chromosomes

var UpDownStreamLen int

var RefIndexStepLen int

var SymbolToEntrezId EntrezIdMap

var RefgeneIndexes RefIndexes

func InitChromArray() {
	refDictFile := path.Join(viper.GetString("db.path"), viper.GetString("db.reference_dict"))
	log.Printf("read %s\n", refDictFile)
	ChromArray = ReadReferenceDictFile(refDictFile)
}

func InitSymbolToEntrezId() {
	symbolToEntrezIdFile := path.Join(viper.GetString("db.path"), viper.GetString("db.gene_symbol_to_entrez_id"))
	log.Printf("read %s\n", symbolToEntrezIdFile)
	SymbolToEntrezId = ReadSymboToIdFile(symbolToEntrezIdFile)
}

func InitRefIndex() {
	refIndexFile := path.Join(viper.GetString("db.path"), viper.GetString("db.refgene_index"))
	log.Printf("read %s\n", refIndexFile)
	RefgeneIndexes = ReadRefIndexFile(refIndexFile)
}

func InitDBParam() {
	UpDownStreamLen = viper.GetInt("param.up_down_stream_length")
	RefIndexStepLen = viper.GetInt("param.refgene_index_step_length")
}

func PrepareDatabase() {
	InitChromArray()
	InitDBParam()
	// reference
	refenceFastaFile := path.Join(viper.GetString("db.path"), viper.GetString("db.reference"))
	log.Printf("read %s\n", refenceFastaFile)
	reference := seq.ReadFastaFile(refenceFastaFile)
	// refgene
	refgeneFile := path.Join(viper.GetString("db.path"), viper.GetString("db.refgene"))
	log.Printf("read %s\n", refgeneFile)
	refgenes := ReadRefgeneFile(refgeneFile)
	// mRNA
	mrnaFastaDir := path.Join(viper.GetString("db.path"), viper.GetString("db.mrna_directory"))
	log.Printf("init and write mRNA to %s\n", mrnaFastaDir)
	if _, err := os.Stat(mrnaFastaDir); os.IsNotExist(err) {
		if err = os.MkdirAll(mrnaFastaDir, os.ModePerm); err != nil {
			log.Panic(err)
		}
	}
	mrnaFasta := NewMrnaFastaMap(refgenes, reference)
	for chrom, fasta := range mrnaFasta {
		mrnaFastaFile := path.Join(mrnaFastaDir, "chr"+chrom+".fa")
		seq.CreateFastaFile(fasta, mrnaFastaFile)
	}
	// Reference Index
	refIndexFile := path.Join(viper.GetString("db.path"), viper.GetString("db.refgene_index"))
	log.Printf("init and write %s\n", refIndexFile)
	CreateRefIndexFile(refgenes, refIndexFile)
}
