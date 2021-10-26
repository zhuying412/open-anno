package snv

import (
	"github.com/spf13/viper"
	"grandanno/db"
	"grandanno/gene"
	"grandanno/seq"
	"log"
	"path"
)

var SplicingDistance int

func InitSnvParam() {
	SplicingDistance = viper.GetInt("param.splicing_distance")
}

func AnnoSnv(inputPath string, outputPath string) {
	// param
	db.InitDBParam()
	InitSnvParam()
	// chrom
	db.InitChromArray()
	// entrez_id
	db.InitSymbolToEntrezId()
	// refgene index
	db.InitRefIndex()
	// refgene
	refgeneFile := path.Join(viper.GetString("db.path"), viper.GetString("db.refgene"))
	log.Printf("read %s\n", refgeneFile)
	allRefgeneMap := gene.ReadRefgeneFile(refgeneFile)
	// avinput
	log.Printf("read %s\n", inputPath)
	allInputs := ReadInputFile(inputPath)
	allResults := make([]Result, 0)
	for _, chrom := range db.ChromArray {
		inputs := allInputs.FilterByChrom(chrom.Name)
		if len(inputs) == 0 {
			continue
		}
		mrnaFastaDir := path.Join(viper.GetString("db.path"), viper.GetString("db.mrna_directory"))
		mrnaFastaFile := path.Join(mrnaFastaDir, "chr"+chrom.Name+".fa")
		log.Printf("read %s\n", mrnaFastaFile)
		mrna := seq.ReadFastaFile(mrnaFastaFile)
		log.Printf("annotate chr%s\n", chrom.Name)
		refgeneMap := allRefgeneMap.FilterByChrom(chrom.Name)
		refIndexes := db.RefgeneIndexes.FilterByChrom(chrom.Name)
		refgeneMap.SetSequence(mrna)
		results := RunAnnotate(inputs, refgeneMap, refIndexes)
		allResults = append(allResults, results...)
	}
	log.Printf("write output %s\n", outputPath)
	CreateAnnotationFile(allResults, outputPath)
}

func init() {
	SplicingDistance = 2
}
