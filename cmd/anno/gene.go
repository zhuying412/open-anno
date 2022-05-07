package anno

import (
	"fmt"
	"log"
	"open-anno/anno/gene"
	"open-anno/pkg/io"
	"open-anno/pkg/io/refgene"
	"open-anno/pkg/seq"
	"os"
	"path"

	"github.com/brentp/faidx"
)

func initGeneBasedData(avinput string, dbPath string, dbName string, builder string) (map[string]io.Variants, refgene.Transcripts, refgene.TransIndexes, map[string]string) {
	// builder
	seq.SetGenome(builder)
	// refgene
	refgeneFile := path.Join(dbPath, builder, dbName, "refgene.txt")
	log.Printf("Read Refgene: %s ...", refgeneFile)
	transcripts, err := refgene.ReadRefgene(refgeneFile)
	if err != nil {
		log.Fatal(err)
	}
	// symbol to id
	ncbiGeneInfoFile := path.Join(dbPath, builder, dbName, "Homo_sapiens.gene_info.gz")
	maneSelectFile := path.Join(dbPath, builder, dbName, "MANE.summary.txt.gz")
	log.Printf("New Gene Symbol to EntrezID from %s, %s", ncbiGeneInfoFile, maneSelectFile)
	symbolToId, err := io.NewGeneSymbolToId(maneSelectFile, ncbiGeneInfoFile, refgeneFile)
	if err != nil {
		log.Fatal(err)
	}
	// index
	indexFile := path.Join(dbPath, builder, dbName, "refgene.idx")
	log.Printf("Read Refgene Index: %s ...", indexFile)
	transIndexes, err := refgene.ReadTransIndexs(indexFile)
	if err != nil {
		log.Fatal(err)
	}
	// snv
	log.Printf("Read avinput: %s ...", avinput)
	snvMap, err := io.ReadVariantMap(avinput)
	if err != nil {
		log.Fatal(err)
	}
	return snvMap, transcripts, transIndexes, symbolToId
}

func initWriter(outfile string) *os.File {
	if outfile == "-" {
		return os.Stdout
	}
	writer, err := os.Create(outfile)
	if err != nil {
		log.Fatal(err)
	}
	return writer
}

func RunAnnoSnvGeneBased(avinput string, dbPath string, dbName string, builder string, outfile string, aashort bool) {
	snvMap, refgenes, refgeneIndexes, symbolToId := initGeneBasedData(avinput, dbPath, dbName, builder)
	// mrna
	mrnaFile := path.Join(dbPath, builder, dbName, "mRNA.fa")
	log.Printf("Read mRNA: %s ...", mrnaFile)
	fai, err := faidx.New(mrnaFile)
	if err != nil {
		log.Fatal(err)
	}
	// anno
	writer := initWriter(outfile)
	fmt.Fprint(writer, "Chr\tStart\tEnd\tRef\tAlt\tGene\tGeneID\tEvent\tRegion\tDetail\n")
	for chrom, snvs := range snvMap {
		log.Printf("Start run annotate chr%s ...", chrom)
		transcripts := refgenes.FilterChrom(chrom, symbolToId, fai, true)
		transIndexes := refgeneIndexes.FilterChrom(chrom)
		if err != nil {
			log.Fatal(err)
		}
		gene.AnnoSnvs(snvs, transcripts, transIndexes, aashort, writer)
	}
}

func RunAnnoCnvGeneBased(avinput string, dbPath string, dbName string, builder string, outfile string) {
	snvMap, refgenes, refgeneIndexes, symbolToId := initGeneBasedData(avinput, dbPath, dbName, builder)
	// anno
	writer := initWriter(outfile)
	fmt.Fprint(writer, "Chr\tStart\tEnd\tRef\tAlt\tRegion\n")
	for chrom, cnvs := range snvMap {
		log.Printf("Filter GeneBased DB by %s ...", chrom)
		transcripts := refgenes.FilterChrom(chrom, symbolToId, &faidx.Faidx{}, false)
		transIndexes := refgeneIndexes.FilterChrom(chrom)
		gene.AnnoCnvs(cnvs, transcripts, transIndexes, writer)
	}
}
