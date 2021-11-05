package snv

import (
	"github.com/spf13/viper"
)

var SplicingDistance int

func InitSnvParam() {
	SplicingDistance = viper.GetInt("param.splicing_distance")
}
