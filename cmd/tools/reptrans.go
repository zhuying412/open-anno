package tools

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"os"
	"path"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type RepTransParam struct {
	RefseqCurated string `validate:"required,pathexists"`
	RefseqHGMD    string `validate:"required,pathexists"`
	RefseqSelect  string `validate:"required,pathexists"`
	RefseqLink    string `validate:"required,pathexists"`
	Output        string `validate:"required"`
}

func (this RepTransParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	validate.Struct(this)
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this RepTransParam) ReadGenePred(gpeFile string, withVer bool) (map[string]map[string]map[string]int, error) {
	data := make(map[string]map[string]map[string]int)
	gpes, err := pkg.ReadGenePred(gpeFile)
	if err != nil {
		return data, err
	}
	for _, gpe := range gpes {
		length := gpe.TxEnd - gpe.TxStart + 1
		chrom, gene, name := gpe.Chrom, gpe.Gene, gpe.Name
		if !withVer {
			name = strings.Split(name, ".")[0]
		}
		if _, ok := data[chrom]; !ok {
			data[chrom] = make(map[string]map[string]int)
		}
		if _, ok := data[chrom][gene]; !ok {
			data[chrom][gene] = make(map[string]int)
		}
		data[chrom][gene][name] = length
	}
	return data, nil
}

func (this RepTransParam) ReadRefseqLink(withVer bool) (map[string]bool, error) {
	data := make(map[string]bool)
	reader, err := pkg.NewIOReader(this.RefseqLink)
	if err != nil {
		return data, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		name := row[0]
		if !withVer {
			name = strings.Split(name, ".")[0]
		}
		data[name] = true
	}
	return data, nil
}

func (this RepTransParam) Run() error {
	curatedData, err := this.ReadGenePred(this.RefseqCurated, true)
	if err != nil {
		return err
	}
	hgmdData, err := this.ReadGenePred(this.RefseqHGMD, false)
	if err != nil {
		return err
	}
	selectData, err := this.ReadGenePred(this.RefseqSelect, false)
	if err != nil {
		return err
	}
	linkData, err := this.ReadRefseqLink(false)
	if err != nil {
		return err
	}
	data := make(map[string]map[string][]string)
	for chrom, subCuratedData1 := range curatedData {
		data[chrom] = make(map[string][]string)
		for gene, subCuratedData2 := range subCuratedData1 {
			data[chrom][gene] = make([]string, 0)
			for trans := range subCuratedData2 {
				name := strings.Split(trans, ".")[0]
				if _, ok := hgmdData[chrom][gene][name]; ok {
					data[chrom][gene] = append(data[chrom][gene], trans)
				}
			}
		}
	}
	for chrom, subCuratedData1 := range curatedData {
		for gene, subCuratedData2 := range subCuratedData1 {
			for trans := range subCuratedData2 {
				if len(data[chrom][gene]) == 0 {
					name := strings.Split(trans, ".")[0]
					if _, ok := selectData[chrom][gene][name]; ok {
						data[chrom][gene] = append(data[chrom][gene], trans)
					}
				}
			}
		}
	}
	for chrom, subCuratedData1 := range curatedData {
		for gene, subCuratedData2 := range subCuratedData1 {
			for trans := range subCuratedData2 {
				if len(data[chrom][gene]) == 0 {
					name := strings.Split(trans, ".")[0]
					if _, ok := linkData[name]; ok {
						data[chrom][gene] = append(data[chrom][gene], trans)
					}
				}
			}
		}
	}
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	for chrom, subData := range data {
		if len(chrom) > 5 {
			continue
		}
		for gene, trans := range subData {
			info := "."
			if len(trans) > 0 {
				info = strings.Join(trans, ",")
			}
			fmt.Fprintf(writer, "%s\t%s\t%s\n", chrom, gene, info)
		}
	}
	return nil
}

func NewRepTransCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "rt",
		Short: "Get Representative Transcript",
		Run: func(cmd *cobra.Command, args []string) {
			var param RepTransParam
			param.RefseqCurated, _ = cmd.Flags().GetString("curated")
			param.RefseqHGMD, _ = cmd.Flags().GetString("hgmd")
			param.RefseqSelect, _ = cmd.Flags().GetString("select")
			param.RefseqLink, _ = cmd.Flags().GetString("link")
			param.Output, _ = cmd.Flags().GetString("output")
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
	cmd.Flags().StringP("curated", "c", "", "Input Refseq Curated File")
	cmd.Flags().StringP("hgmd", "g", "", "Input Refseq HGMD File")
	cmd.Flags().StringP("select", "s", "", "Input Refseq Select File")
	cmd.Flags().StringP("link", "l", "", "Input Refseq Link File")
	cmd.Flags().StringP("output", "o", "", "Output Merged Annotation File")
	return cmd
}
