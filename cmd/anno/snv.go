package anno

import (
	"log"
	"open-anno/anno/db"
	"open-anno/anno/gene/snv"
	"open-anno/pkg/schema"
	"path"
	"strings"

	"github.com/spf13/pflag"
)

type AnnoSnvGBParam struct {
	AnnoParam
	AAshort  bool
	Exon     bool
	WithAggs bool
}

func (this AnnoSnvGBParam) Mrna() string {
	return path.Join(this.DBpath, this.Builder, this.DBname+"_mRNA.fa")
}

func (this AnnoSnvGBParam) GeneData() (snv.GeneData, error) {
	return snv.NewGeneData(this.DBFile(), this.DBIndex(), this.GeneInfo(), this.Mrna())
}

func (this *AnnoSnvGBParam) Bind(flagSet *pflag.FlagSet) {
	this.AnnoParam.Bind(flagSet)
	this.WithAggs, _ = flagSet.GetBool("with_aggs")
	this.AAshort, _ = flagSet.GetBool("aashort")
	this.Exon, _ = flagSet.GetBool("exon")
}

func (this AnnoSnvGBParam) Valid() error {
	schema.IS_EXON_REGION = this.Exon
	return this.AnnoParam.Valid()
}

func (this AnnoSnvGBParam) Run() error {
	geneData, err := this.GeneData()
	if err != nil {
		return err
	}
	if this.WithAggs {
		outdir, outnames := path.Dir(this.Output), strings.Split(path.Base(this.Output), ".")
		if len(outnames) == 1 {
			outnames = append(outnames, "raw")
		} else {
			outnames = append(outnames[0:len(outnames)-1], "raw", outnames[len(outnames)-1])
		}
		tempOutput := path.Join(outdir, strings.Join(outnames, "."))
		err = snv.AnnoSnvs(this.Input, tempOutput, this.DBname, geneData, this.AAshort)
		if err != nil {
			return err
		}
		log.Println("Run Aggregating ...")
		return snv.AggsTransAnno(tempOutput, this.DBname, this.Output)
	}
	return snv.AnnoSnvs(this.Input, this.Output, this.DBname, geneData, this.AAshort)
}

type AnnoSnvFBParam struct {
	AnnoParam
}

func (this AnnoSnvFBParam) Run() error {
	return db.AnnoFilterBased(this.Input, this.DBFile(), this.Output)
}

type AnnoSnvRBParam struct {
	AnnoParam
	Overlap float64 `validate:"omitempty,min=0,max=1"`
}

func (this *AnnoSnvRBParam) Bind(flagSet *pflag.FlagSet) {
	this.AnnoParam.Bind(flagSet)
	this.Overlap, _ = flagSet.GetFloat64("overlap")
}

func (this AnnoSnvRBParam) Run() error {
	return db.AnnoRegionBased(this.Input, this.DBFile(), this.Output, this.Overlap)
}
