package main

import (
	"flag"
	"fmt"
	"grandanno/config"
	"grandanno/prepare"
	"path"
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
	//referenceFile := path.Join(param.DatabasePath, config.MyConfig.Database.Reference)
	refgeneFile := path.Join(param.DatabasePath, config.MyConfig.Database.Refgene)
	ensMtFile := path.Join(param.DatabasePath, config.MyConfig.Database.EnsMt)
	//mrnaFile := path.Join(param.DatabasePath, config.MyConfig.Database.Mrna)
	refidxFile := path.Join(param.DatabasePath, config.MyConfig.Database.Refidx)
	upDownSteamLen := config.MyConfig.Param.UpDownStream
	refidxStepLen := config.MyConfig.Param.RefidxStep
	refgeneDict := make(prepare.RefgeneDict)
	refgeneDict.Read(refgeneFile, upDownSteamLen)
	refgeneDict.Read(ensMtFile, upDownSteamLen)
	//refgeneDict.Write(referenceFile, mrnaFile)
	refidxDict := make(prepare.RefidxDict)
	refidxDict.Init(refidxStepLen)
	refidxDict.Write(refidxFile, refgeneDict)
}

func main() {
	//var param Parameter
	//param.Init()
	//switch {
	//case param.Help:
	//	flag.Usage()
	//case param.IsPrepare:
	//	param.runPrepare()
	//}
	b := make([]int, 10, -1)
	fmt.Print(b)
}
