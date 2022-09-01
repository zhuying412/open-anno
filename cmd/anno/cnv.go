package anno

import (
	"open-anno/anno/gene/cnv"
)

type AnnoCnvGBParam struct {
	AnnoParam
}

func (this AnnoCnvGBParam) GeneData() (cnv.GeneData, error) {
	return cnv.NewGeneData(this.DBFile(), this.DBIndex(), this.GeneInfo())
}

func (this AnnoCnvGBParam) Run() error {
	geneData, err := this.GeneData()
	if err != nil {
		return err
	}
	return cnv.AnnoCnvs(this.Input, this.Output(), this.DBname, geneData)
}

type AnnoCnvRBParam struct {
	AnnoSnvRBParam
}
