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

func initGeneBasedData(avinput, dbPath, dbName, builder string) (map[string]io.Variants, refgene.Transcripts, refgene.TransIndexes, map[string]string, error) {
	var snvMap map[string]io.Variants
	var transcripts refgene.Transcripts
	var transIndexes refgene.TransIndexes
	var symbolToId map[string]string
	var err error
	// builder
	seq.SetGenome(builder)
	// refgene
	refgeneFile := path.Join(dbPath, builder, dbName, "refgene.txt")
	log.Printf("Read Refgene: %s ...", refgeneFile)
	transcripts, err = refgene.ReadRefgene(refgeneFile)
	if err != nil {
		return snvMap, transcripts, transIndexes, symbolToId, err
	}
	// symbol to id
	ncbiGeneInfoFile := path.Join(dbPath, builder, "Homo_sapiens.gene_info.gz")
	gene2RefseqFile := path.Join(dbPath, builder, "Homo_sapiens.gene2refseq.gz")
	log.Printf("New Gene Symbol to EntrezID from %s, %s", ncbiGeneInfoFile, gene2RefseqFile)
	symbolToId, err = io.NewGeneSymbolToId(gene2RefseqFile, ncbiGeneInfoFile, refgeneFile)
	if err != nil {
		return snvMap, transcripts, transIndexes, symbolToId, err
	}
	// index
	indexFile := path.Join(dbPath, builder, dbName, "refgene.idx")
	log.Printf("Read Refgene Index: %s ...", indexFile)
	transIndexes, err = refgene.ReadTransIndexs(indexFile)
	if err != nil {
		return snvMap, transcripts, transIndexes, symbolToId, err
	}
	// snv
	log.Printf("Read avinput: %s ...", avinput)
	snvMap, err = io.ReadVariantMap(avinput)
	if err != nil {
		return snvMap, transcripts, transIndexes, symbolToId, err
	}
	return snvMap, transcripts, transIndexes, symbolToId, err
}

func initWriter(outfile string) (io.WriteCloser, error) {
	if outfile == "-" {
		return os.Stdout, nil
	}
	return io.NewIoWriter(outfile)
}

func RunAnnoSnvGeneBased(avinput string, dbPath string, dbName string, builder string, outfile string, aashort bool, errChan chan error) {
	snvMap, refgenes, refgeneIndexes, symbolToId, err := initGeneBasedData(avinput, dbPath, dbName, builder)
	if err != nil {
		errChan <- err
		return
	}
	// mrna
	mrnaFile := path.Join(dbPath, builder, dbName, "mRNA.fa")
	log.Printf("Read mRNA: %s ...", mrnaFile)
	fai, err := faidx.New(mrnaFile)
	if err != nil {
		errChan <- err
		return
	}
	// anno
	writer, err := initWriter(outfile)
	if err != nil {
		errChan <- err
		return
	}
	fmt.Fprint(writer, "Chr\tStart\tEnd\tRef\tAlt\tGene\tGeneID\tEvent\tRegion\tDetail\n")
	for chrom, snvs := range snvMap {
		log.Printf("Start run annotate chr%s ...", chrom)
		transcripts := refgenes.FilterChrom(chrom, symbolToId, fai, true)
		transIndexes := refgeneIndexes.FilterChrom(chrom)
		if err != nil {
			errChan <- err
			return
		}
		gene.AnnoSnvs(snvs, transcripts, transIndexes, aashort, writer)
	}
	errChan <- err
}

func RunAnnoCnvGeneBased(avinput string, dbPath string, dbName string, builder string, outfile string, errChan chan error) {
	snvMap, refgenes, refgeneIndexes, symbolToId, err := initGeneBasedData(avinput, dbPath, dbName, builder)
	if err != nil {
		errChan <- err
		return
	}
	// anno
	writer, err := initWriter(outfile)
	if err != nil {
		errChan <- err
		return
	}
	fmt.Fprint(writer, "Chr\tStart\tEnd\tRef\tAlt\tRegion\n")
	for chrom, cnvs := range snvMap {
		log.Printf("Filter GeneBased DB by %s ...", chrom)
		transcripts := refgenes.FilterChrom(chrom, symbolToId, &faidx.Faidx{}, false)
		transIndexes := refgeneIndexes.FilterChrom(chrom)
		gene.AnnoCnvs(cnvs, transcripts, transIndexes, writer)
	}
	errChan <- err
}
