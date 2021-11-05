package db

import (
	"github.com/spf13/viper"
	"log"
	"path"
)

var ChromArray Chromosomes

var UpDownStreamLen int

var RefIndexStepLen int

var DBIndexStepLen int

var SymbolToEntrezId EntrezIdMap

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

func InitDBParam() {
	UpDownStreamLen = viper.GetInt("param.up_down_stream_length")
	RefIndexStepLen = viper.GetInt("param.refgene_index_step_length")
}

func GetRefIndexes() RefIndexes {
	refIndexFile := path.Join(viper.GetString("db.path"), viper.GetString("db.refgene_index"))
	log.Printf("read %s\n", refIndexFile)
	return ReadRefIndexFile(refIndexFile)
}

func GetMrnaFile(chrom Chromosome) string {
	mrnaFastaDir := path.Join(viper.GetString("db.path"), viper.GetString("db.mrna_directory"))
	return path.Join(mrnaFastaDir, "chr"+chrom.Name+".fa")
}

func init() {
	DBIndexStepLen = 100000
}
