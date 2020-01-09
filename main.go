package main

import (
	"flag"
	"grandanno/prepare"
)

type Parameter struct {
	IsPrepare      bool
	DatabasePath   string
	SplicingLength int
	Help           bool
}

func (param *Parameter) Init() {
	flag.BoolVar(&param.IsPrepare, "p", false, "Is Prepare Database")
	flag.StringVar(&param.DatabasePath, "d", "./database", "Directory Path of Database")
	flag.IntVar(&param.SplicingLength, "s", 15, "Length of Splicing Region")
	flag.BoolVar(&param.Help, "h", false, "Help")
	flag.Parse()
}

func (param Parameter) runPrepare() {
	var refgeneDict prepare.RefgeneDict
	refgeneDict.Read(param.DatabasePath)
}

func main() {
	var param Parameter
	param.Init()
	switch {
	case param.Help:
		flag.Usage()
	case param.IsPrepare:
		param.runPrepare()
	}
}
