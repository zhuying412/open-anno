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
	reader, err := pkg.NewIOReader(this.GenePred)
	if err != nil {
		return err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		trans, err := pkg.NewTranscript(scanner.Text())
		if err != nil {
			return err
		}
		if _, ok := repTransData[trans.PK()]; ok {
			if this.CDS {
				cdss := make([][]int, 0)
				for i := 0; i < trans.ExonCount; i++ {
					start, end := trans.ExonStarts[i], trans.ExonEnds[i]
					if trans.CdsStart <= end && start <= trans.CdsEnd {
						if start < trans.CdsStart && trans.CdsStart < end {
							start = trans.CdsStart
						}
						if start < trans.CdsEnd && trans.CdsEnd < end {
							end = trans.CdsEnd
						}
						cdss = append(cdss, []int{start, end})
					}
				}
				cdsCount := len(cdss)
				for i := 0; i < cdsCount; i++ {
					cds := fmt.Sprintf("CDS%d", i+1)
					if trans.Strand == "-" {
						cds = fmt.Sprintf("CDS%d", cdsCount-i)
					}
					fmt.Fprintf(writer, "%s\t%d\t%d\t%s:%s:%s\n", trans.Chrom, cdss[i][0]-1, cdss[i][1], trans.Gene, trans.Name, cds)
				}
			} else {
				for i := 0; i < trans.ExonCount; i++ {
					exon := fmt.Sprintf("Exon%d", i+1)
					if trans.Strand == "-" {
						exon = fmt.Sprintf("Exon%d", trans.ExonCount-i)
					}
					fmt.Fprintf(writer, "%s\t%d\t%d\t%s:%s:%s\n", trans.Chrom, trans.ExonStarts[i]-1, trans.ExonEnds[i], trans.Gene, trans.Name, exon)
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
