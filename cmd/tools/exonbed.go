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

type ExonBedParam struct {
	GenePred string `validate:"required,pathexists"`
	RepTrans string `validate:"required,pathexists"`
	CDS      bool
	Output   string `validate:"required"`
}

func (this ExonBedParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this ExonBedParam) ReadRepTrans() (map[string]bool, error) {
	data := make(map[string]bool)
	reader, err := pkg.NewIOReader(this.RepTrans)
	if err != nil {
		return data, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		key := strings.Join(row[0:3], ":")
		data[key] = true
	}
	return data, nil
}

func (this ExonBedParam) Run() error {
	repTransData, err := this.ReadRepTrans()
	if err != nil {
		return err
	}
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	gpes, err := pkg.ReadGenePred(this.GenePred)
	for _, gpe := range gpes {
		key := fmt.Sprintf("%s:%s:%s", gpe.Chrom, gpe.Gene, gpe.Name)
		if _, ok := repTransData[key]; ok {
			if this.CDS {
				cdss := make([][]int, 0)
				for i := 0; i < gpe.ExonCount; i++ {
					start, end := gpe.ExonStarts[i], gpe.ExonEnds[i]
					if gpe.CdsStart <= end && start <= gpe.CdsEnd {
						if start < gpe.CdsStart && gpe.CdsStart < end {
							start = gpe.CdsStart
						}
						if start < gpe.CdsEnd && gpe.CdsEnd < end {
							end = gpe.CdsEnd
						}
						cdss = append(cdss, []int{start, end})
					}
				}
				cdsCount := len(cdss)
				for i := 0; i < cdsCount; i++ {
					cds := fmt.Sprintf("CDS%d", i+1)
					if gpe.Strand == "-" {
						cds = fmt.Sprintf("CDS%d", cdsCount-i)
					}
					fmt.Fprintf(writer, "%s\t%d\t%d\t%s:%s:%s\n", gpe.Chrom, cdss[i][0]-1, cdss[i][1], gpe.Gene, gpe.Name, cds)
				}
			} else {
				for i := 0; i < gpe.ExonCount; i++ {
					exon := fmt.Sprintf("Exon%d", i+1)
					if gpe.Strand == "-" {
						exon = fmt.Sprintf("Exon%d", gpe.ExonCount-i)
					}
					fmt.Fprintf(writer, "%s\t%d\t%d\t%s:%s:%s\n", gpe.Chrom, gpe.ExonStarts[i]-1, gpe.ExonEnds[i], gpe.Gene, gpe.Name, exon)
				}
			}
		}
	}
	return nil
}

func NewExonBedCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "eb",
		Short: "Get Exon BED from GenePred",
		Run: func(cmd *cobra.Command, args []string) {
			var param ExonBedParam
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.RepTrans, _ = cmd.Flags().GetString("transcript")
			param.Output, _ = cmd.Flags().GetString("output")
			param.CDS, _ = cmd.Flags().GetBool("cds")
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
	cmd.Flags().StringP("genepred", "g", "", "Input GenePed File, such as Refseq Curated from UCSC")
	cmd.Flags().StringP("transcript", "t", "", "Input Representative Transcript File")
	cmd.Flags().BoolP("cds", "c", false, "Parameter, Output CDS BED instead of Exon")
	cmd.Flags().StringP("output", "o", "", "Output BED File")
	return cmd
}
