package gene

import (
	"github.com/spf13/viper"
	"log"
	"path"
)

var SymbolToEntrezId EntrezIdMap

func Init() {
	symbolToEntrezIdFile := path.Join(viper.GetString("db.path"), viper.GetString("db.gene_entrez"))
	log.Printf("read %s\n", symbolToEntrezIdFile)
	SymbolToEntrezId = ReadSymboToIdFile(symbolToEntrezIdFile)
}
