package index

import (
	"OpenAnno/db/chromosome"
	"log"
	"os"
	"path"
)

func Generate(refgeneFile string, transcriptDir string, upDownStreamLen int, indexStepLen int) {
	err := os.MkdirAll(transcriptDir, os.ModePerm)
	if err != nil {
		log.Panic(err)
	}
	log.Printf("read %s", refgeneFile)
	transcripts := ReadRefgeneFile(refgeneFile, upDownStreamLen)
	log.Print("init transcript indexes")
	indexes := InitTranscriptIndexes(indexStepLen)
	for _, chrom := range chromosome.ChromList {
		outfile := path.Join(transcriptDir, "chr"+chrom.Name+".idx.json")
		_, err := os.Stat(outfile)
		if os.IsExist(err) {
			continue
		}
		log.Printf("write %s", outfile)
		subTranscripts := transcripts.FilterByChrom(chrom.Name)
		subIndexes := indexes.FilterByChrom(chrom.Name)
		if subTranscripts.Len() > 0 && subIndexes.Len() > 0 {
			CreateTranscriptIndexJSON(subIndexes, subTranscripts, outfile)
		}
	}
}
