package gene

import (
	"github.com/spf13/viper"
	"log"
	"path"
)

var SymbolToEntrezId EntrezIdMap

func Init() {
	if len(SymbolToEntrezId) == 0 {
		symbolToEntrezIdFile := path.Join(viper.GetString("db.path"), viper.GetString("db.gene_entrez"))
		log.Printf("read %s", symbolToEntrezIdFile)
		SymbolToEntrezId = ReadSymboToIdFile(symbolToEntrezIdFile)
	}
}
