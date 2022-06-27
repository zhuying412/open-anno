package anno

import (
	"errors"
	"fmt"
	"open-anno/anno/db"
	"open-anno/anno/gene"
	"open-anno/anno/gene/cnv"
	"open-anno/anno/gene/snv"
	"open-anno/pkg/io/refgene"
	"open-anno/pkg/seq"
	"os"
	"path"

	"github.com/go-playground/validator/v10"
)

func CheckPathExists(fl validator.FieldLevel) bool {
	path := fl.Field().String()
	_, err := os.Stat(path)
	return os.IsExist(err)
}

type AnnoParam struct {
	Input     string `validate:"required,pathexists"`
	DBType    string `validate:"required,oneof=g f r"`
	DBpath    string `validate:"required,pathexists"`
	DBname    string `validate:"required,pathexists"`
	Builder   string `validate:"required"`
	OutPrefix string `validate:"required"`
	AAshort   bool
	Exon      bool
	Overlap   float64 `validate:"omitempty,min=0,max=1"`
}

func (this AnnoParam) Output() string {
	return fmt.Sprintf("%s.%s.anno.txt", this.OutPrefix, this.DBname)
}

func (this AnnoParam) DBFile() string {
	return path.Join(this.DBpath, this.Builder, this.DBname+".txt")
}

func (this AnnoParam) DBIndex() string {
	return path.Join(this.DBpath, this.Builder, this.DBname+".txt.idx")
}

func (this AnnoParam) Mrna() string {
	return path.Join(this.DBpath, this.Builder, this.DBname+"_mRNA.fa")
}

func (this AnnoParam) NcbiGeneInfo() string {
	return path.Join(this.DBpath, this.Builder, "Homo_sapiens.gene_info.gz")
}

func (this AnnoParam) Gene2Refseq() string {
	return path.Join(this.DBpath, this.Builder, "Homo_sapiens.gene2refseq.gz")
}

func (this AnnoParam) GeneData() (gene.GeneData, error) {
	if this.DBType == "g" {
		return gene.NewGeneData(this.DBFile(), this.DBFile(), this.NcbiGeneInfo(), this.Gene2Refseq(), this.Mrna())
	}
	return gene.GeneData{}, errors.New("Non GeneBased has no GeneData")
}

func (this AnnoParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	seq.SetGenome(this.Builder)
	refgene.IS_EXON_REGION = this.Exon
	return nil
}

func (this AnnoParam) RunAnno(varType string, errChan chan error) {
	err := this.Valid()
	if err != nil {
		errChan <- err
		return
	}
	if this.DBType == "g" {
		geneData, err := this.GeneData()
		if err != nil {
			errChan <- err
			return
		}
		if varType == "snv" {
			errChan <- snv.AnnoSnvs(this.Input, this.Output(), this.DBname, geneData, this.AAshort)
		}
		if varType == "cnv" {
			errChan <- cnv.AnnoCnvs(this.Input, this.Output(), this.DBname, geneData)
		}
	}
	if this.DBType == "f" {
		errChan <- db.AnnoFilterBased(this.Input, this.DBFile(), this.Output())
	}
	if this.DBType == "r" {
		errChan <- db.AnnoRegionBased(this.Input, this.DBFile(), this.Output(), this.Overlap)
	}
}
