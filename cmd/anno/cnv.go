package anno

import (
	"fmt"
	"log"
	"open-anno/anno"
	"open-anno/anno/db"
	"open-anno/anno/gene"
	"open-anno/pkg"
	"os"
	"path"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type AnnoCnvParam struct {
	Input         string   `validate:"required,pathexists"`
	GenePred      string   `validate:"required,pathexists"`
	GenePredIndex string   `validate:"required,pathexists"`
	Gene          string   `validate:"required,pathexists"`
	Output        string   `validate:"required"`
	GBName        string   `validate:"required"`
	RegionBaseds  []string `validate:"pathsexists"`
	Overlap       float64  `validate:"required"`
	Clean         bool
}

func (this *AnnoCnvParam) Valid() error {
	this.GenePredIndex = this.GenePred + ".idx"
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	validate.RegisterValidation("pathsexists", pkg.CheckPathsExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return os.MkdirAll(this.Outdir(), 0666)
}

func (this AnnoCnvParam) Outdir() string {
	return path.Dir(this.Output)
}

func (this *AnnoCnvParam) RunAnnoGeneBase(cnvs anno.Variants, annoChan chan anno.AnnoInfos, headerChan chan string, errChan chan error) {
	// 读取GenePred输入文件
	log.Printf("Read GenePred: %s ...", this.GenePred)
	gpes, err := pkg.ReadGenePred(this.GenePred)
	if err != nil {
		errChan <- err
		return
	}
	// 读取GenePred Index输入文件
	log.Printf("Read GenePred Index: %s ...", this.GenePredIndex)
	transIndexes, err := pkg.ReadTransIndexes(this.GenePredIndex)
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
	annoInfos, err := gene.AnnoCnvs(cnvs, gpes, transIndexes, geneSymbolToID)
	annoChan <- annoInfos
	headerChan <- "Region"
	errChan <- err
}

func (this *AnnoCnvParam) RunAnnoRegionBase(cnvs anno.Variants, database string, annoChan chan anno.AnnoInfos, headerChan chan string, errChan chan error) {
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
	annoInfos, err := db.AnnoRegionBased(cnvs, regVars, scanner.Name, this.Overlap)
	annoChan <- annoInfos
	headerChan <- scanner.Name
	errChan <- err
}

func (this AnnoCnvParam) Run() error {
	// 读取变异输入文件
	log.Printf("Read AnnoInput: %s ...", this.Input)
	cnvs, err := anno.ReadBED(this.Input)
	if err != nil {
		return err
	}
	// 构造 error channel
	annoChan := make(chan anno.AnnoInfos, len(this.RegionBaseds)+1)
	errChan := make(chan error, len(this.RegionBaseds)+1)
	headerChan := make(chan string, len(this.RegionBaseds)+1)
	// GeneBased 注释
	go this.RunAnnoGeneBase(cnvs, annoChan, headerChan, errChan)
	// RegionBased 注释
	for _, database := range this.RegionBaseds {
		go this.RunAnnoRegionBase(cnvs, database, annoChan, headerChan, errChan)
	}
	// 错误处理
	annoInfos := make(anno.AnnoInfos)
	for i := 0; i < len(this.RegionBaseds)+1; i++ {
		err := <-errChan
		if err != nil {
			return err
		}
		for pk, infos := range <-annoChan {
			if items, ok := annoInfos[pk]; ok {
				annoInfos[pk] = append(items, infos...)
			} else {
				annoInfos[pk] = infos
			}
		}
	}
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprint(writer, "Chrom\tStart\tRef\tAlt\tAnnotation\n")
	for _, cnv := range cnvs {
		annoVar := cnv.AnnoVariant()
		annoTexts := make([]string, 0)
		for _, annoInfo := range annoInfos[annoVar.PK()] {
			annoTexts = append(annoTexts, fmt.Sprintf("%s=%s", annoInfo.Key, annoInfo.Value))
		}
		fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s%s\n", annoVar.Chrom, annoVar.Start, annoVar.End, annoVar.Ref, annoVar.Alt, strings.Join(annoTexts, ";"))
	}
	return nil
}

func NewAnnoCnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "cnv",
		Short: "Annotate CNV",
		Run: func(cmd *cobra.Command, args []string) {
			var param AnnoCnvParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.Gene, _ = cmd.Flags().GetString("gene")
			param.GBName, _ = cmd.Flags().GetString("gbname")
			param.Output, _ = cmd.Flags().GetString("output")
			param.RegionBaseds, _ = cmd.Flags().GetStringArray("regionbaseds")
			param.Overlap, _ = cmd.Flags().GetFloat64("overlap")
			param.Clean, _ = cmd.Flags().GetBool("clean")
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
	cmd.Flags().StringP("gene", "g", "", "Input Gene Symbol To ID File")
	cmd.Flags().StringP("gbname", "n", "", "Parameter Database Name")
	cmd.Flags().StringP("output", "o", "", "AnnoOutput File")
	cmd.Flags().StringArrayP("regionbaseds", "r", []string{}, "Input RegionBased Database File")
	cmd.Flags().Float64P("overlap", "l", 0.7, "Parameter Database Name")
	cmd.Flags().BoolP("clean", "c", false, "Clean Temporary File")
	return cmd
}
