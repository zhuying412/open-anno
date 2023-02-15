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
	ClinvarHGVS   string `validate:"required,pathexists"`
	Output        string `validate:"required"`
}

func (this RepTransParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
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

func (this RepTransParam) ReadClinvarHGVS() (map[string][]string, error) {
	data := make(map[string][]string)
	reader, err := pkg.NewIOReader(this.ClinvarHGVS)
	if err != nil {
		return data, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		if !strings.HasPrefix(row[0], "#") &&
			row[0] != "-" && (row[5] == "GRCh38" || row[5] == "na") &&
			(strings.HasPrefix(row[6], "NM_") || strings.HasPrefix(row[6], "NR_")) &&
			row[10] == "Yes" && row[11] != "No" {
			gene := row[0]
			trans := strings.Split(row[6], ".")[0]
			if _, ok := data[gene]; !ok {
				data[gene] = make([]string, 0)
			}
			if pkg.FindArr(data[gene], trans) == -1 {
				data[gene] = append(data[gene], trans)
			}
		}
	}
	return data, nil
}

// func (this RepTransParam) ReadClinvarHGVS() (map[string]string, error) {
// 	data := make(map[string]string)
// 	reader, err := pkg.NewIOReader(this.ClinvarHGVS)
// 	if err != nil {
// 		return data, err
// 	}
// 	defer reader.Close()
// 	tmpData := make(map[string]map[string]int)
// 	scanner := pkg.NewIOScanner(reader)
// 	for scanner.Scan() {
// 		row := strings.Split(scanner.Text(), "\t")
// 		if !strings.HasPrefix(row[0], "#") && row[0] != "-" && (row[5] == "GRCh38" || row[5] == "na") && (strings.HasPrefix(row[6], "NM_") || strings.HasPrefix(row[6], "NR_")) && row[10] == "Yes" {
// 			gene := row[0]
// 			trans := strings.Split(row[6], ".")[0]
// 			if _, ok := tmpData[gene]; !ok {
// 				tmpData[gene] = make(map[string]int)
// 			}
// 			if _, ok := tmpData[gene][trans]; !ok {
// 				tmpData[gene][trans] = 0
// 			}
// 			tmpData[gene][trans]++
// 		}
// 	}
// 	for gene, transcripts := range tmpData {
// 		maxCount := 0
// 		maxCountTrans := ""
// 		for trans, count := range transcripts {
// 			if maxCount < count {
// 				maxCount = count
// 				maxCountTrans = trans
// 			}
// 		}
// 		data[gene] = maxCountTrans
// 	}
// 	return data, nil
// }

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
	clinvarData, err := this.ReadClinvarHGVS()
	if err != nil {
		return err
	}
	data := make(map[string]map[string]string)
	maxLenData := make(map[string]map[string]string)
	for chrom, genes := range curatedData {
		data[chrom] = make(map[string]string)
		maxLenData[chrom] = make(map[string]string)
		for gene, transcripts := range genes {
			data[chrom][gene] = "."
			maxLen := 0
			maxLenTrans := ""
			for trans, length := range transcripts {
				if maxLen < length {
					maxLen = length
					maxLenTrans = trans
				}
			}
			maxLenData[chrom][gene] = maxLenTrans
		}
	}

	for chrom, genes := range curatedData {
		for gene, transcripts := range genes {
			if data[chrom][gene] == "." {
				for trans := range transcripts {
					name := strings.Split(trans, ".")[0]
					if _, ok := hgmdData[chrom][gene][name]; ok {
						data[chrom][gene] = trans + "\tHGMD"
						break
					}
				}
			}
		}
	}
	for chrom, genes := range curatedData {
		for gene, transcripts := range genes {
			if data[chrom][gene] == "." {
				for trans := range transcripts {
					name := strings.Split(trans, ".")[0]
					if _, ok := selectData[chrom][gene][name]; ok {
						data[chrom][gene] = trans + "\tSelect"
						break
					}
				}
			}
		}
	}
	for chrom, genes := range curatedData {
		for gene, transcripts := range genes {
			if data[chrom][gene] == "." {
				for trans := range transcripts {
					if transArr, ok := clinvarData[gene]; ok {
						name := strings.Split(trans, ".")[0]
						if pkg.FindArr(transArr, name) != -1 {
							data[chrom][gene] = trans + "\tClinVar"
							break
						}
					}
					// if trans == clinvarData[gene] {
					// 	data[chrom][gene] = trans + "\tClinVar"
					// 	break
					// }
				}
			}
		}
	}
	for chrom, genes := range curatedData {
		for gene := range genes {
			if data[chrom][gene] == "." {
				data[chrom][gene] = maxLenData[chrom][gene] + "\tMaxLength"
			}
		}
	}
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	for chrom, genes := range data {
		for gene, trans := range genes {
			fmt.Fprintf(writer, "%s\t%s\t%s\n", chrom, gene, trans)
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
			param.ClinvarHGVS, _ = cmd.Flags().GetString("clinvar")
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
	cmd.Flags().StringP("curated", "c", "", "Input Refseq Curated File from UCSC, name: ncbiRefSeqCurated.txt.gz")
	cmd.Flags().StringP("hgmd", "g", "", "Input Refseq HGMD File  from UCSC, name: ncbiRefSeqHgmd.txt.gz")
	cmd.Flags().StringP("select", "s", "", "Input Refseq Select File from UCSC, name: ncbiRefSeqSelect.txt.gz")
	cmd.Flags().StringP("clinvar", "l", "", "Input ClinVar HGVS File from NCBI, name: hgvs4variation.txt.gz")
	cmd.Flags().StringP("output", "o", "", "Output Representative Transcript File")
	return cmd
}
