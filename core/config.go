package core

import (
	"errors"
	"gopkg.in/yaml.v2"
	"io/ioutil"
)

type YamlData interface {
	GetYaml(yamlDir string) (yamlFile string)
}

type ConfFile struct {
	Reference string `yaml:"reference"`
	NcbiGene  string `yaml:"ncbi_gene"`
	Refgene   string `yaml:"refgene"`
	EnsMt     string `yaml:"ens_mt"`
	Cds       string `yaml:"cds"`
	Exon      string `yaml:"exon"`
	Mrna      string `yaml:"mrna"`
	Refidx    string `yaml:"refidx"`
}

type ConfParam struct {
	UpDownStream int `yaml:"up_down_stream"`
	RefidxStep   int `yaml:"refidx_step"`
	SplicingLen  int `yaml:"splicing_len"`
}

type ConfChrom []struct {
	Name   string `yaml:"name"`
	Length int    `yaml:"length"`
}

func (c ConfChrom) GetNames() (chromList []string) {
	for _, chrom := range c {
		chromList = append(chromList, chrom.Name)
	}
	return
}

func (c ConfChrom) GetByName(name string) (order int, length int) {
	for index, chrom := range c {
		if order, length = index+1, chrom.Length; chrom.Name == name {
			return order, length
		}
	}
	panic(errors.New("Not Found: " + name))
}

type Config struct {
	File  ConfFile  `yaml:"file"`
	Param ConfParam `yaml:"param"`
	Chrom ConfChrom `yaml:"chrom"`
}

func (c *Config) Read(yamlFile string) {
	buffer, _ := ioutil.ReadFile(yamlFile)
	if err := yaml.Unmarshal(buffer, c); err != nil {
		panic(err)
	}
}

var Conf Config
