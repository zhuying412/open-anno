package anno

import (
	"fmt"
	"log"
	"open-anno/anno/db"
	"open-anno/pkg/io"
	"open-anno/pkg/seq"
	"os"
	"path"
)

func RunAnnoFilterBased(avinput string, dbPath string, dbName string, builder string, outfile string, errChan chan error) {
	// builder
	seq.SetGenome(builder)
	// snvs
	log.Printf("Read avinput: %s ...", avinput)
	snvs, err := io.ReadVariantMap(avinput)
	if err != nil {
		errChan <- err
		return
	}
	// writer
	var writer io.WriteCloser = os.Stdout
	if outfile != "-" {
		writer, err = io.NewIoWriter(outfile)
		if err != nil {
			errChan <- err
			return
		}

	}
	// anno
	headerWrited := true
	for chrom, subSnvs := range snvs {
		log.Printf("Start run annotate %s chr%s ...", dbName, chrom)
		dbFile := path.Join(dbPath, builder, dbName, fmt.Sprintf("chr%s.txt", chrom))
		if _, err = os.Stat(dbFile); os.IsNotExist(err) {
			continue
		}
		db.AnnoFilterBased(subSnvs, dbFile, headerWrited, writer)
		headerWrited = false
	}
	errChan <- err
}

func RunAnnoRegionBased(avinput string, dbPath string, dbName string, builder string, overlap float64, outfile string, errChan chan error) {
	// builder
	seq.SetGenome(builder)
	// snvs
	log.Printf("Read avinput: %s ...", avinput)
	snvs, err := io.ReadVariantMap(avinput)
	if err != nil {
		errChan <- err
		return
	}
	// writer
	var writer io.WriteCloser = os.Stdout
	if outfile != "-" {
		writer, err = io.NewIoWriter(outfile)
		if err != nil {
			errChan <- err
			return
		}

	}
	// anno
	headerWrited := true
	for chrom, subSnvs := range snvs {
		dbFile := path.Join(dbPath, builder, dbName, fmt.Sprintf("chr%s.txt", chrom))
		if _, err = os.Stat(dbFile); os.IsNotExist(err) {
			continue
		}
		db.AnnoRegionBased(subSnvs, dbFile, overlap, headerWrited, writer)
		headerWrited = false
	}
	errChan <- err
}
