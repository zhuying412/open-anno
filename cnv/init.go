package cnv

import (
	"github.com/spf13/viper"
	"grandanno/db"
	"grandanno/gene"
	"log"
	"path"
)

func AnnoCnv(inputPath string, outputPath string) {
	// param
	db.InitDBParam()
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
	// input
	log.Printf("read %s\n", inputPath)
	allInputs := ReadInputFile(inputPath)
	allResults := make([]Result, 0)
	for _, chrom := range db.ChromArray {
		inputs := allInputs.FilterByChrom(chrom.Name)
		if len(inputs) == 0 {
			continue
		}
		log.Printf("annotate chr%s\n", chrom.Name)
		refgeneMap := allRefgeneMap.FilterByChrom(chrom.Name)
		refIndexes := db.RefgeneIndexes.FilterByChrom(chrom.Name)
		results := RunAnnotate(inputs, refgeneMap, refIndexes)
		allResults = append(allResults, results...)
	}
	log.Printf("write output %s\n", outputPath)
	CreateAnnotationFile(allResults, outputPath)
}
