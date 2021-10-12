package cnv

import (
	"grandanno/db"
	"grandanno/gene"
	"grandanno/seq"
	"log"
)

func AnnoCnv(databasePath string, genomeVersion string, avinputPath string, outputPath string) {
	db.DB = db.NewDatabase(databasePath, genomeVersion)
	db.ChromArray = db.NewChromosomesFromReferenceDict(db.DB.ReferenceDictFile())
	// mRNA
	log.Printf("start read %s\n", db.DB.MrnaFastaFile())
	mrna := seq.NewFasta(db.DB.MrnaFastaFile())
	log.Printf("start read %s\n", db.DB.RefIndexFile())
	refIndexes := db.NewRefIndexes(db.DB.RefIndexFile())
	log.Printf("start read %s\n", db.DB.RefgeneFile())
	refgeneMap := gene.NewRefgeneMap(db.DB.RefgeneFile())
	refgeneMap.SetSequence(mrna)
	log.Printf("start read %s\n", avinputPath)
	snvs := NewCnvs(avinputPath)
	log.Print("run annotate")
	annos := RunAnnotate(snvs, refgeneMap, refIndexes)
	log.Printf("write outputPath %s\n", avinputPath)
	WriteAnnotations(annos, outputPath)
}
