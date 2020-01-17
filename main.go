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
	VarType        string
	DatabasePath   string
	SplicingLength int
	Help           bool
}

func (param *Parameter) Init() {
	flag.BoolVar(&param.IsPrepare, "p", false, "Is Prepare Database")
	flag.StringVar(&param.Input, "i", "", "Input File")
	flag.StringVar(&param.Ouput, "o", "", "Output File")
	flag.StringVar(&param.VarType, "t", "gatk_snv", "The Source of variants:[gatk_snv, xhmm_cnv]")
	flag.StringVar(&param.DatabasePath, "d", "./database", "Directory Path of Database")
	flag.IntVar(&param.SplicingLength, "s", 15, "Length of Splicing Region")
	flag.BoolVar(&param.Help, "h", false, "Help")
	flag.Parse()
}

func (param Parameter) runPrepare() {
	// File Path
	referenceFile := path.Join(param.DatabasePath, config.MyConfig.Database.Reference)
	refgeneFile := path.Join(param.DatabasePath, config.MyConfig.Database.Refgene)
	ensMtFile := path.Join(param.DatabasePath, config.MyConfig.Database.EnsMt)
	mrnaFile := path.Join(param.DatabasePath, config.MyConfig.Database.Mrna)
	refidxFile := path.Join(param.DatabasePath, config.MyConfig.Database.Refidx)
	// Param
	upDownSteamLen := config.MyConfig.Param.UpDownStream
	refidxStepLen := config.MyConfig.Param.RefidxStep
	// mRNA
	reference := make(core.Fasta)
	reference.Read(referenceFile)
	refgeneDict := make(prepare.RefgeneDict)
	refgeneDict.Read(refgeneFile, upDownSteamLen)
	refgeneDict.Read(ensMtFile, upDownSteamLen)
	refgeneDict.Write(reference, mrnaFile)
	// refidx
	refidxDict := make(prepare.RefidxDict)
	refidxDict.Init(refidxStepLen)
	refidxDict.Write(refidxFile, refgeneDict)
}

func (param Parameter) runAnnoGatkSnv() {
	// File Path
	ncbiGeneInfoFile := path.Join(param.DatabasePath, config.MyConfig.Database.NcbiGene)
	mrnaFile := path.Join(param.DatabasePath, config.MyConfig.Database.Mrna)
	refgeneFile := path.Join(param.DatabasePath, config.MyConfig.Database.Refgene)
	ensMtFile := path.Join(param.DatabasePath, config.MyConfig.Database.EnsMt)
	refidxFile := path.Join(param.DatabasePath, config.MyConfig.Database.Refidx)
	// Param
	upDownSteamLen := config.MyConfig.Param.UpDownStream
	splicingLen := config.MyConfig.Param.SplicingLen
	// NCBI Gene Info
	ncbiGene := core.NcbiGene{}
	ncbiGene.Read(ncbiGeneInfoFile)
	// Refgene
	mrna := make(core.Fasta)
	mrna.Read(mrnaFile)
	refgeneDict := make(core.RefgeneDict)
	refgeneDict.Read(refgeneFile, ncbiGene)
	refgeneDict.Read(ensMtFile, ncbiGene)
	refgeneDict.SetSequence(mrna, upDownSteamLen)
	refidxs := make(core.Refidxs, 0)
	refidxs.Read(refidxFile)
	// VCF
	vcf := snv.Vcf{File: param.Input}
	snvs := vcf.ReadAll()
	// Annotate
	reults := make(snv.Results, 0)
	reults.RunAnno(snvs, refgeneDict, refidxs, splicingLen)
	reults.Write(param.Ouput)

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
	case param.VarType == "gatk_snv":
		param.runAnnoGatkSnv()
	default:
		param.runAnnoGatkSnv()
	}
}

func main() {
	a, b := "ATGC", "ATCC"
	var i int
	for i = 0; i < len(a); i++ {
		if a[i] != b[i] {
			break
		}
	}
	fmt.Println()
}
