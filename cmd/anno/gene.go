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
	"strings"

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
	refgeneFile := path.Join(dbPath, builder, dbName, "GenePred.txt")
	log.Printf("Read GenePred: %s ...", refgeneFile)
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
	indexFile := path.Join(dbPath, builder, dbName, "GenePred.idx")
	log.Printf("Read GenePred Index: %s ...", indexFile)
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
	// writer
	var writer io.WriteCloser = os.Stdout
	if !strings.HasPrefix(outfile, "-") {
		writer, err = io.NewIoWriter(outfile)
		if err != nil {
			errChan <- err
			return
		}

	}
	fmt.Fprintf(writer,
		"Chr\tStart\tEnd\tRef\tAlt\t%s.Gene\t%s.GeneID\t%s.Event\t%s.Region\t%s.Detail\n",
		dbName, dbName, dbName, dbName, dbName,
	)
	// anno
	for chrom, snvs := range snvMap {
		log.Printf("Start run annotate %s chr%s ...", dbName, chrom)
		transcripts, err := refgenes.FilterChrom(chrom, symbolToId, fai, true)
		if err != nil {
			errChan <- err
			return
		}
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
	cnvMap, refgenes, refgeneIndexes, symbolToId, err := initGeneBasedData(avinput, dbPath, dbName, builder)
	if err != nil {
		errChan <- err
		return
	}
	// writer
	var writer io.WriteCloser = os.Stdout
	if !strings.HasPrefix(outfile, "-") {
		writer, err = io.NewIoWriter(outfile)
		if err != nil {
			errChan <- err
			return
		}

	}
	fmt.Fprintf(writer, "Chr\tStart\tEnd\tRef\tAlt\t%s.Region\n", dbName)
	// anno
	for chrom, cnvs := range cnvMap {
		log.Printf("Filter GeneBased DB by %s ...", chrom)
		transcripts, err := refgenes.FilterChrom(chrom, symbolToId, &faidx.Faidx{}, false)
		if err != nil {
			errChan <- err
			return
		}
		transIndexes := refgeneIndexes.FilterChrom(chrom)
		gene.AnnoCnvs(cnvs, transcripts, transIndexes, writer)
	}
	errChan <- err
}
