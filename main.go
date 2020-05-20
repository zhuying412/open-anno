package main

import (
	"flag"
	"fmt"
	"grandanno/cnv"
	"grandanno/core"
	"grandanno/prepare"
	"grandanno/snv"
	"os"
	"path"
)

type Parameter struct {
	Input          string
	Ouput          string
	Config         string
	IsPrepare      bool
	VarType        string
	DatabasePath   string
	SplicingLength int
	Help           bool
}

func (param Parameter) runPrepare() {
	// mRNA
	fmt.Printf("start read %s\n", core.Conf.File.Reference)
	reference := make(core.Fasta)
	reference.Read(core.Conf.File.Reference)
	fmt.Printf("start read %s and %s\n", core.Conf.File.Refgene, core.Conf.File.EnsMt)
	refgeneDict := make(prepare.RefgeneDict)
	refgeneDict.Read(core.Conf.File.Refgene, core.Conf.Param.UpDownStream)
	refgeneDict.Read(core.Conf.File.EnsMt, core.Conf.Param.UpDownStream)
	fmt.Printf("start write %s\n", core.Conf.File.Mrna)
	refgeneDict.Write(reference, core.Conf.File.Mrna)
	// Reference Index
	fmt.Printf("start init and write %s\n", core.Conf.File.Refidx)
	refidxDict := make(prepare.RefidxDict)
	refidxDict.Init(core.Conf.Param.RefidxStep)
	refidxDict.Write(core.Conf.File.Refidx, refgeneDict)
}

func (param Parameter) runAnnoGatkSnv() {
	// NCBI Gene Info
	fmt.Printf("start read %s\n", core.Conf.File.NcbiGene)
	ncbiGene := core.NcbiGene{}
	ncbiGene.Read(core.Conf.File.NcbiGene)
	// Refgene
	fmt.Printf("start read %s\n", core.Conf.File.Mrna)
	mrna := make(core.Fasta)
	mrna.Read(core.Conf.File.Mrna)
	fmt.Printf("start read %s and %s\n", core.Conf.File.Refgene, core.Conf.File.EnsMt)
	refgeneDict := make(core.RefgeneDict)
	refgeneDict.Read(core.Conf.File.Refgene, ncbiGene)
	refgeneDict.Read(core.Conf.File.EnsMt, ncbiGene)
	refgeneDict.SetSequence(mrna)
	refgeneDict.SetUpDownStream(core.Conf.Param.UpDownStream)
	refidxs := make(core.Refidxs, 0)
	refidxs.Read(core.Conf.File.Refidx)
	// VCF
	fmt.Printf("start read %s\n", param.Input)
	snvs := make(snv.Snvs, 0)
	snvs.ReadGatkVcfFile(param.Input)
	// Annotate
	fmt.Printf("start annoate and write out %s\n", param.Ouput)
	results := make(snv.Results, 0)
	results.RunAnno(snvs, refgeneDict, refidxs, core.Conf.Param.SplicingLen)
	results.Write(param.Ouput)
}

func (param Parameter) runAnnoXhmmCnv() {
	// NCBI Gene Info
	ncbiGene := core.NcbiGene{}
	ncbiGene.Read(core.Conf.File.NcbiGene)
	// Refgene
	refgeneDict := make(core.RefgeneDict)
	refgeneDict.Read(core.Conf.File.Refgene, ncbiGene)
	refgeneDict.Read(core.Conf.File.EnsMt, ncbiGene)
	refgeneDict.SetUpDownStream(core.Conf.Param.UpDownStream)
	refidxs := make(core.Refidxs, 0)
	refidxs.Read(core.Conf.File.Refidx)
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

var param Parameter

func init() {
	flag.BoolVar(&param.IsPrepare, "p", false, "Is Prepare Database")
	flag.StringVar(&param.Input, "i", "", "Input File")
	flag.StringVar(&param.Ouput, "o", "", "Output File")
	flag.StringVar(&param.Config, "c", "", "Config YAML File")
	flag.StringVar(&param.VarType, "t", "gatk_snv", "The Source of variants:[gatk_snv, xhmm_cnv]")
	flag.StringVar(&param.DatabasePath, "d", "./database", "Directory Path of Database")
	flag.IntVar(&param.SplicingLength, "s", 15, "Length of Splicing Region")
	flag.BoolVar(&param.Help, "h", false, "Help")
	flag.Parse()
	if !param.IsPrepare && (len(param.Input) == 0 || len(param.Ouput) == 0 || len(param.Config) == 0) {
		flag.Usage()
		os.Exit(-1)
	}
	core.Conf.Read(param.Config)
}

func main() {
	switch {
	case param.Help:
		flag.Usage()
	case param.IsPrepare:
		param.runPrepare()
	case param.VarType == "gatk_snv":
		param.runAnnoGatkSnv()
	case param.VarType == "xhmm_snv":
		param.runAnnoXhmmCnv()
	default:
		flag.Usage()
	}
}
