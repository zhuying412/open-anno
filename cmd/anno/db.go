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

func RunAnnoFilterBased(avinput string, dbPath string, dbName string, builder string, outfile string) {
	// builder
	seq.SetGenome(builder)
	// snvs
	log.Printf("Read avinput: %s ...", avinput)
	snvs, err := io.ReadVariantMap(avinput)
	if err != nil {
		log.Fatal(err)
	}
	// anno
	var writer *os.File
	if outfile == "-" {
		writer = os.Stdout
	} else {
		writer, err = os.Create(outfile)
		if err != nil {
			log.Fatal(err)
		}
	}
	headerWrited := true
	for chrom, subSnvs := range snvs {
		dbFile := path.Join(dbPath, builder, dbName, fmt.Sprintf("chr%s.txt", chrom))
		if _, err = os.Stat(dbFile); os.IsNotExist(err) {
			continue
		}
		db.AnnoFilterBased(subSnvs, dbFile, headerWrited, writer)
		headerWrited = false
	}
}

func RunAnnoRegionBased(avinput string, dbPath string, dbName string, builder string, overlap float64, outfile string) {
	// builder
	seq.SetGenome(builder)
	// snvs
	log.Printf("Read avinput: %s ...", avinput)
	snvs, err := io.ReadVariantMap(avinput)
	if err != nil {
		log.Fatal(err)
	}
	// anno
	var writer *os.File
	if outfile == "-" {
		writer = os.Stdout
	} else {
		writer, err = os.Create(outfile)
		if err != nil {
			log.Fatal(err)
		}
	}
	headerWrited := true
	for chrom, subSnvs := range snvs {
		dbFile := path.Join(dbPath, builder, dbName, fmt.Sprintf("chr%s.txt", chrom))
		if _, err = os.Stat(dbFile); os.IsNotExist(err) {
			continue
		}
		db.AnnoRegionBased(subSnvs, dbFile, overlap, headerWrited, writer)
		headerWrited = false
	}
}
