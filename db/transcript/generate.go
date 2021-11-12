package transcript

import (
	"OpenAnno/db/chromosome"
	"OpenAnno/seq"
	"log"
	"os"
	"path"
)

func Generate(refgeneFile string, referenceFile string, transcriptDir string, upDownStreamLen int) {
	//Init()
	err := os.MkdirAll(transcriptDir, os.ModePerm)
	if err != nil {
		log.Panic(err)
	}
	log.Printf("read %s\n", referenceFile)
	reference := seq.ReadFastaFile(referenceFile)
	log.Printf("read %s\n", refgeneFile)
	transcripts := ReadRefgeneFile(refgeneFile, upDownStreamLen)
	for _, chrom := range chromosome.ChromList {
		outfile := path.Join(transcriptDir, "chr"+chrom.Name+".json")
		_, err := os.Stat(outfile)
		if os.IsExist(err) {
			continue
		}
		log.Printf("write %s\n", outfile)
		subTranscripts := transcripts.FilterByChrom(chrom.Name)
		chromSeq := reference[chrom.Name]
		if subTranscripts.Len() > 0 {
			CreateTranscriptJSON(subTranscripts, chromSeq, outfile)
		}
	}
}
