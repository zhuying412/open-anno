package main

import (
	"flag"
	"fmt"
	"grandanno/cnv"
	"grandanno/config"
	"grandanno/core"
	"grandanno/prepare"
	"grandanno/snv"
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
	fmt.Printf("start read %s\n", referenceFile)
	reference := make(core.Fasta)
	reference.Read(referenceFile)
	fmt.Printf("start read %s and %s\n", refgeneFile, ensMtFile)
	refgeneDict := make(prepare.RefgeneDict)
	refgeneDict.Read(refgeneFile, upDownSteamLen)
	refgeneDict.Read(ensMtFile, upDownSteamLen)
	fmt.Printf("start write %s\n", mrnaFile)
	refgeneDict.Write(reference, mrnaFile)
	// refidx
	fmt.Printf("start init and write %s\n", refidxFile)
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
	fmt.Printf("start read %s\n", ncbiGeneInfoFile)
	ncbiGene := core.NcbiGene{}
	ncbiGene.Read(ncbiGeneInfoFile)
	// Refgene
	fmt.Printf("start read %s\n", mrnaFile)
	mrna := make(core.Fasta)
	mrna.Read(mrnaFile)
	fmt.Printf("start read %s and %s\n", refgeneFile, ensMtFile)
	refgeneDict := make(core.RefgeneDict)
	refgeneDict.Read(refgeneFile, ncbiGene)
	refgeneDict.Read(ensMtFile, ncbiGene)
	refgeneDict.SetSequence(mrna)
	refgeneDict.SetUpDownStream(upDownSteamLen)
	refidxs := make(core.Refidxs, 0)
	refidxs.Read(refidxFile)
	// VCF
	fmt.Printf("start read %s\n", param.Input)
	snvs := make(snv.Snvs, 0)
	snvs.ReadGatkVcfFile(param.Input)
	// Annotate
	fmt.Printf("start annoate and write out %s\n", param.Ouput)
	results := make(snv.Results, 0)
	results.RunAnno(snvs, refgeneDict, refidxs, splicingLen)
	results.Write(param.Ouput)
}

func (param Parameter) runAnnoXhmmCnv() {
	// File Path
	ncbiGeneInfoFile := path.Join(param.DatabasePath, config.MyConfig.Database.NcbiGene)
	refgeneFile := path.Join(param.DatabasePath, config.MyConfig.Database.Refgene)
	ensMtFile := path.Join(param.DatabasePath, config.MyConfig.Database.EnsMt)
	refidxFile := path.Join(param.DatabasePath, config.MyConfig.Database.Refidx)
	// Param
	upDownSteamLen := config.MyConfig.Param.UpDownStream
	// NCBI Gene Info
	ncbiGene := core.NcbiGene{}
	ncbiGene.Read(ncbiGeneInfoFile)
	// Refgene
	refgeneDict := make(core.RefgeneDict)
	refgeneDict.Read(refgeneFile, ncbiGene)
	refgeneDict.Read(ensMtFile, ncbiGene)
	refgeneDict.SetUpDownStream(upDownSteamLen)
	refidxs := make(core.Refidxs, 0)
	refidxs.Read(refidxFile)
	// VCF
	cnvDict := make(cnv.XhmmCnvDict, 0)
	cnvDict.ReadXhmmVcfFile(param.Input)
	// Annotate
	for sample, cnvs := range cnvDict {
		outJsonFile := path.Join(param.Ouput + "." + sample + ".json")
		results := make(cnv.Results, 0)
		results.RunAnno(cnvs, refgeneDict, refidxs)
		results.Write(outJsonFile)
	}
}

//func main() {
//	var param Parameter
//	param.Init()
//	if !param.IsPrepare && (len(param.Input) == 0 || len(param.Ouput) == 0) {
//		flag.Usage()
//		os.Exit(-1)
//	}
//	switch {
//	case param.Help:
//		flag.Usage()
//	case param.IsPrepare:
//		param.runPrepare()
//	case param.VarType == "gatk_snv":
//		param.runAnnoGatkSnv()
//	case param.VarType == "xhmm_snv":
//		param.runAnnoXhmmCnv()
//	default:
//		flag.Usage()
//	}
//}

type SS string

func (s *SS) test() {
	*s = "a"
}
func main() {
	var s SS
	s.test()
	fmt.Println(s)
}
