package config

import (
	"gopkg.in/yaml.v2"
	"io/ioutil"
)

type Database struct {
	Reference string
	NcbiGene  string
	Refgene   string
	EnsMt     string
	Cds       string
	Exon      string
	Mrna      string
	Refidx    string
}

type Param struct {
	UpDownStream int
}
type Config struct {
	Database Database
	Param    Param
}

func (config *Config) Init(yamlFile string) {
	buffer, _ := ioutil.ReadFile(yamlFile)
	if err := yaml.Unmarshal(buffer, &config); err != nil {
		panic(err.Error())
	}
}

func readConfig() Config {
	yamlFile := "./config/config.yaml"
	config := Config{}
	config.Init(yamlFile)
	return config
}

var (
	MyConfig = readConfig()
)
