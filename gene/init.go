package gene

import (
	"github.com/spf13/viper"
	"log"
	"path"
)

func GetRefgeneMap() RefgeneMap {
	refgeneFile := path.Join(viper.GetString("db.path"), viper.GetString("db.refgene"))
	log.Printf("read %s\n", refgeneFile)
	return ReadRefgeneFile(refgeneFile)
}
