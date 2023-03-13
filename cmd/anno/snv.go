package anno

import (
	"log"
	"open-anno/anno"
	"open-anno/anno/db"
	"open-anno/anno/gene"
	"open-anno/pkg"
	"os"
	"path"

	"github.com/brentp/faidx"
	"github.com/brentp/vcfgo"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type AnnoSnvParam struct {
	Input              string `validate:"required,pathexists"`
	GenePred           string `validate:"required,pathexists"`
	GenePredIndex      string `validate:"required,pathexists"`
	Genome             string `validate:"required,pathexists"`
	GenomeIndex        string `validate:"required,pathexists"`
	Gene               string `validate:"required,pathexists"`
	Output             string `validate:"required"`
	AAshort            bool
	Exon               bool
	FilterBaseds       []string `validate:"pathsexists"`
	FilterBasedIndexes []string `validate:"pathsexists"`
	RegionBaseds       []string `validate:"pathsexists"`
	Overlap            float64  `validate:"required"`
}

func (this *AnnoSnvParam) Valid() error {
	this.GenePredIndex = this.GenePred + ".idx"
	this.GenomeIndex = this.Genome + ".fai"
	for _, db := range this.FilterBaseds {
		this.FilterBasedIndexes = append(this.FilterBasedIndexes, db+".idx")
	}
	pkg.IS_EXON_REGION = this.Exon
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	validate.RegisterValidation("pathsexists", pkg.CheckPathsExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return os.MkdirAll(this.Outdir(), 0666)
}

func (this AnnoSnvParam) Outdir() string {
	return path.Dir(this.Output)
}

func (this *AnnoSnvParam) RunAnnoGeneBase(snvs anno.Variants, annoChan chan anno.AnnoInfos, vcfHeaderInfoChan chan map[string]*vcfgo.Info, errChan chan error) {
	// 读取GenePred输入文件
	log.Printf("Read GenePred: %s ...", this.GenePred)
	gpes, err := pkg.ReadGenePred(this.GenePred)
	if err != nil {
		errChan <- err
		return
	}
	// 读取Transcript Index输入文件
	log.Printf("Read GenePred Index: %s ...", this.GenePredIndex)
	transIndexes, err := pkg.ReadTransIndexes(this.GenePredIndex)
	if err != nil {
		errChan <- err
		return
	}
	// 读取Genome输入文件
	log.Printf("Read Reference Fasta: %s ...", this.Genome)
	genome, err := faidx.New(this.Genome)
	if err != nil {
		errChan <- err
		return
	}
	// 读取gene输入文件
	log.Printf("Read Gene: %s ...", this.Gene)
	geneSymbolToID, err := anno.ReadGene(this.Gene)
	if err != nil {
		errChan <- err
		return
	}
	annoInfos, err := gene.AnnoSnvs(snvs, gpes, transIndexes, genome, geneSymbolToID, this.AAshort)
	vcfHeaderInfoChan <- map[string]*vcfgo.Info{
		"GENE": {
			Id:          "GENE",
			Description: "Gene Symbol",
			Number:      ".",
			Type:        "String",
		},
		"GENE_ID": {
			Id:          "GENE_ID",
			Description: "Gene Entrez ID",
			Number:      ".",
			Type:        "String",
		},
		"REGION": {
			Id:          "REGION",
			Description: "Region in gene, eg: exonic, intronic, UTR3, UTR5",
			Number:      ".",
			Type:        "String",
		},
		"EVENT": {
			Id:          "EVENT",
			Description: "Variant Event, eg: missense, nonsense, splicing",
			Number:      ".",
			Type:        "String",
		},
		"DETAIL": {
			Id:          "DETAIL",
			Description: "Gene detail, FORMAT=Gene:Transcript:Exon:NA_CHANGE:AA_CHANGE",
			Number:      ".",
		},
	}
	annoChan <- annoInfos
	errChan <- err
}

func (this *AnnoSnvParam) RunAnnoFilterBase(snvs anno.Variants, database, databaseIndex string, annoChan chan anno.AnnoInfos, vcfHeaderInfoChan chan map[string]*vcfgo.Info, errChan chan error) {
	reader, err := db.NewFilterVarReader(database)
	if err != nil {
		errChan <- err
		return
	}
	defer reader.Close()
	idxs, binSize, err := db.ReadFilterVarIdx(databaseIndex)
	if err != nil {
		errChan <- err
		return
	}

	annoInfos, err := db.AnnoFilterBased(snvs, reader, idxs, binSize)
	annoChan <- annoInfos
	vcfHeaderInfoChan <- reader.VCFHeaderInfo()
	errChan <- err
}

func (this *AnnoSnvParam) RunAnnoRegionBase(snvs anno.Variants, database string, annoChan chan anno.AnnoInfos, vcfHeaderInfoChan chan map[string]*vcfgo.Info, errChan chan error) {
	reader, err := pkg.NewIOReader(database)
	if err != nil {
		errChan <- err
		return
	}
	defer reader.Close()
	scanner := db.NewRegionVarScanner(reader)
	regVars, err := scanner.ReadAll()
	if err != nil {
		errChan <- err
		return
	}
	annoInfos, err := db.AnnoRegionBased(snvs, regVars, scanner.Name, this.Overlap)
	annoChan <- annoInfos
	vcfHeaderInfoChan <- scanner.VCFHeaderInfo()
	errChan <- err
}

func (this AnnoSnvParam) Run() error {
	// 读取变异输入文件
	log.Printf("Read AnnoInput: %s ...", this.Input)
	snvs, vcfHeader, err := anno.ReadVCF(this.Input)
	if err != nil {
		return err
	}
	// 构造 error channel
	errChan := make(chan error, len(this.FilterBaseds)+len(this.RegionBaseds)+1)
	annoInfosChan := make(chan anno.AnnoInfos, len(this.FilterBaseds)+len(this.RegionBaseds)+1)
	vcfHeaderInfoChan := make(chan map[string]*vcfgo.Info, len(this.FilterBaseds)+len(this.RegionBaseds)+1)
	go this.RunAnnoGeneBase(snvs, annoInfosChan, vcfHeaderInfoChan, errChan)
	// FilterBased 注释
	for i, database := range this.FilterBaseds {
		go this.RunAnnoFilterBase(snvs, database, this.FilterBasedIndexes[i], annoInfosChan, vcfHeaderInfoChan, errChan)
	}
	// RegionBased 注释
	for _, database := range this.RegionBaseds {
		go this.RunAnnoRegionBase(snvs, database, annoInfosChan, vcfHeaderInfoChan, errChan)
	}
	// 错误处理
	annoInfos := make(anno.AnnoInfos)
	for i := 0; i < len(this.FilterBaseds)+len(this.RegionBaseds)+1; i++ {
		err := <-errChan
		if err != nil {
			return err
		}
		for pk, infos := range <-annoInfosChan {
			if items, ok := annoInfos[pk]; ok {
				annoInfos[pk] = append(items, infos...)
			} else {
				annoInfos[pk] = infos
			}
		}
		for key, info := range <-vcfHeaderInfoChan {
			vcfHeader.Infos[key] = info
		}
	}
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()

	vcfWriter, err := vcfgo.NewWriter(writer, vcfHeader)
	for _, snv := range snvs {
		vcfVariant := snv.(anno.SNV).Variant
		for _, annoInfo := range annoInfos[snv.AnnoVariant().PK()] {
			if annoInfo.Value != "" && annoInfo.Value != "." {
				err = vcfVariant.Info().Set(annoInfo.Key, annoInfo.Value)
				if err != nil {
					return err
				}
			}
		}
		vcfWriter.WriteVariant(&vcfVariant)
	}
	return nil
}

func NewAnnoSnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "snv",
		Short: "Annotate  SNV",
		Run: func(cmd *cobra.Command, args []string) {
			var param AnnoSnvParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.Genome, _ = cmd.Flags().GetString("genome")
			param.Gene, _ = cmd.Flags().GetString("gene")
			param.Output, _ = cmd.Flags().GetString("output")
			param.AAshort, _ = cmd.Flags().GetBool("aashort")
			param.Exon, _ = cmd.Flags().GetBool("exon")
			param.FilterBaseds, _ = cmd.Flags().GetStringArray("filterbaseds")
			param.RegionBaseds, _ = cmd.Flags().GetStringArray("regionbaseds")
			param.Overlap, _ = cmd.Flags().GetFloat64("overlap")
			err := param.Valid()
			if err != nil {
				cmd.Help()
				log.Fatal(err)
			}
			err = param.Run()
			if err != nil {
				log.Fatal(err)
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "AnnoInput File")
	cmd.Flags().StringP("genepred", "d", "", "Input GenePred File")
	cmd.Flags().StringP("genome", "G", "", "Input Genome Fasta File")
	cmd.Flags().StringP("gene", "g", "", "Input Gene Symbol To ID File")
	cmd.Flags().StringP("output", "o", "", "AnnoOutput File")
	cmd.Flags().BoolP("aashort", "a", false, "Parameter Is AA Short")
	cmd.Flags().BoolP("exon", "e", false, "Parameter Is Exon")
	cmd.Flags().StringArrayP("filterbaseds", "f", []string{}, "Input FilterBased Database File")
	cmd.Flags().StringArrayP("regionbaseds", "r", []string{}, "Input RegionBased Database File")
	cmd.Flags().Float64P("overlap", "l", 0.7, "Parameter Database Name")
	return cmd
}
