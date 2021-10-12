package db

import (
	"log"
	"os"
	"path"
)

type Database struct {
	Path          string
	GenomeVersion string
}

func (d Database) IsExists(dbPath string) bool {
	_, err := os.Stat(dbPath)
	return !os.IsNotExist(err)
}

func (d Database) ReferenceFastaFile() string {
	return path.Join(d.Path, d.GenomeVersion+".fa")
}

func (d Database) ReferenceDictFile() string {
	return path.Join(d.Path, d.GenomeVersion+".dict")
}

func (d Database) GeneSymbolToIDFile() string {
	return path.Join(d.Path, "GeneSymbolToID")
}

func (d Database) RefgeneFile() string {
	return path.Join(d.Path, d.GenomeVersion+"_refGene_ensMT.txt")
}

func (d Database) RefIndexFile() string {
	return path.Join(d.Path, d.GenomeVersion+"_refgene_ensMT.idx")
}

func (d Database) CdsBedFile() string {
	return path.Join(d.Path, d.GenomeVersion+"_cds.bed")
}

func (d Database) ExonBedFile() string {
	return path.Join(d.Path, d.GenomeVersion+"_exon.bed")
}

func (d Database) MrnaFastaFile() string {
	return path.Join(d.Path, d.GenomeVersion+"_mRNA.fa")

}

func (d Database) Check() {
	//log.Printf("Reference Fasta: %s", d.ReferenceFasta())
	for _, dbPath := range []string{
		d.ReferenceDictFile(), d.GeneSymbolToIDFile(), d.RefgeneFile(),
		d.RefIndexFile(), d.CdsBedFile(), d.ExonBedFile(), d.MrnaFastaFile(),
	} {
		if !d.IsExists(dbPath) {
			log.Printf("%s Exists", dbPath)
		} else {
			log.Printf("%s Not found", dbPath)
		}
	}
}

func NewDatabase(databaePath string, genomeVersion string) Database {
	return Database{Path: databaePath, GenomeVersion: genomeVersion}
}
