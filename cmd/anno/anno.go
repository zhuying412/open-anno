package anno

import (
	"errors"
	"fmt"
	"log"
	"open-anno/anno/db"
	"open-anno/anno/gene"
	"open-anno/anno/gene/cnv"
	"open-anno/anno/gene/snv"
	"open-anno/cmd/pre"
	"open-anno/pkg/io/refgene"
	"open-anno/pkg/seq"
	"path"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type AnnoParam struct {
	Input     string `validate:"required,pathexists"`
	DBType    string `validate:"required,oneof=g f r"`
	DBpath    string `validate:"required,pathexists"`
	DBname    string `validate:"required"`
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
		return gene.NewGeneData(this.DBFile(), this.DBIndex(), this.NcbiGeneInfo(), this.Gene2Refseq(), this.Mrna())
	}
	return gene.GeneData{}, errors.New("Non GeneBased has no GeneData")
}

func (this AnnoParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pre.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	seq.SetGenome(this.Builder)
	refgene.IS_EXON_REGION = this.Exon
	return nil
}

func (this AnnoParam) RunAnno(varType string, errChan chan error) {
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

func NewAnnoParams(cmd *cobra.Command) ([]AnnoParam, error) {
	input, _ := cmd.Flags().GetString("avinput")
	outprefix, _ := cmd.Flags().GetString("outprefix")
	dbpath, _ := cmd.Flags().GetString("dbpath")
	builder, _ := cmd.Flags().GetString("builder")
	aashort, _ := cmd.Flags().GetBool("aashort")
	exon, _ := cmd.Flags().GetBool("exon")
	overlap, _ := cmd.Flags().GetFloat64("overlap")
	dbtypes, _ := cmd.Flags().GetString("dbtypes")
	dbnames, _ := cmd.Flags().GetString("dbnames")
	dbNames := strings.Split(dbnames, ",")
	dbTypes := strings.Split(dbtypes, ",")
	annoParams := make([]AnnoParam, len(dbNames))
	for i := 0; i < len(dbNames); i++ {
		dbname := strings.TrimSpace(dbNames[i])
		dbtype := strings.TrimSpace(dbTypes[i])
		param := AnnoParam{
			Input:     input,
			OutPrefix: outprefix,
			DBpath:    dbpath,
			Builder:   builder,
			AAshort:   aashort,
			Exon:      exon,
			Overlap:   overlap,
			DBname:    dbname,
			DBType:    dbtype,
		}
		err := param.Valid()
		if err != nil {
			cmd.Help()
			return annoParams, err
		}
		annoParams[i] = param
	}
	return annoParams, nil
}

func NewAnnoCmd(varType string) *cobra.Command {
	varType = strings.ToLower(varType)
	if varType != "snv" && varType != "cnv" {
		log.Fatalln("only 'snv' or 'cnv'")
	}
	cmd := &cobra.Command{
		Use:   varType,
		Short: fmt.Sprintf("Annotate for %s", strings.ToUpper(varType)),
		Run: func(cmd *cobra.Command, args []string) {
			annoParams, err := NewAnnoParams(cmd)
			if err != nil {
				log.Fatal(err)
			}
			errChan := make(chan error, len(annoParams))
			for _, annoParam := range annoParams {
				go annoParam.RunAnno(varType, errChan)
			}
			for i := 0; i < len(annoParams); i++ {
				err := <-errChan
				if err != nil {
					log.Fatal(err)
				}
			}
		},
	}
	cmd.Flags().StringP("avinput", "i", "", "Annotated Variants Input File")
	cmd.Flags().StringP("outprefix", "o", "", "Output Prefix")
	cmd.Flags().StringP("dbpath", "d", "", "Database Directory")
	cmd.Flags().StringP("dbnames", "n", "", "Database Names")
	cmd.Flags().StringP("dbtypes", "t", "", "Database Types")
	cmd.Flags().StringP("builder", "b", "hg19", "Database Builder")
	if varType == "snv" {
		cmd.Flags().BoolP("aashort", "s", false, "Database Builder")
		cmd.Flags().BoolP("exon", "e", false, "Output ExonOrder Instead of TypeOrder")
	}
	if varType == "cnv" {
		cmd.Flags().Float64P("overlap", "p", 0.75, "CNV Overlap Threshold")
	}
	return cmd
}
