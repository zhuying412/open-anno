package main

import (
	"flag"
	"fmt"
	"grandanno/config"
	"grandanno/core"
	"grandanno/prepare"
	"grandanno/snv"
	"os"
	"path"
)

type Parameter struct {
	Input string
	Ouput string
	//Type string
	IsPrepare      bool
	DatabasePath   string
	SplicingLength int
	Help           bool
}

func (param *Parameter) Init() {
	flag.BoolVar(&param.IsPrepare, "p", false, "Is Prepare Database")
	flag.StringVar(&param.Input, "i", "", "Input File")
	flag.StringVar(&param.Ouput, "o", "", "Output File")
	//flag.StringVar(&param.Type, "t", "gatk", "Input Type:[gatk, xhmm]")
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

func (param Parameter) runAnnoGatk() {
	// NCBI Gene Info
	var ncbiGene core.NcbiGene
	ncbiGeneInfoFile := path.Join(param.DatabasePath, config.MyConfig.Database.NcbiGene)
	ncbiGene.Read(ncbiGeneInfoFile)
	// Refgene
	// VCF
	// Annotate
	vcf := snv.Vcf{File: param.Input}
	snvs := vcf.ReadAll()
	for _, s := range snvs {
		fmt.Println(s)
	}
}

func main1() {
	var param Parameter
	param.Init()
	if len(param.Input) == 0 || len(param.Ouput) == 0 {
		flag.Usage()
		os.Exit(-1)
	}
	switch {
	case param.Help:
		flag.Usage()
	case param.IsPrepare:
		param.runPrepare()
	default:
		param.runAnnoGatk()
	}
}
func main() {
	main1()
	//a,b := 3,5
	//
	//fmt.Print(float32(a)/float32(b))
}
