package clinvar

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"os"
	"path"
	"strings"

	"github.com/brentp/vcfgo"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type PreClinvarGeneParam struct {
	Input  string `validate:"required,pathexists"`
	Output string `validate:"required"`
}

func (this PreClinvarGeneParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this PreClinvarGeneParam) Run() error {
	reader, err := pkg.NewIOReader(this.Input)
	if err != nil {
		return err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return err
	}
	defer vcfReader.Close()
	counts := make(map[string][]int)
	for row := vcfReader.Read(); row != nil; row = vcfReader.Read() {
		row.Chromosome = "chr" + row.Chromosome
		imc, err := row.Info().Get("MC")
		if err != nil {
			continue
		}
		igeneinfo, err := row.Info().Get("GENEINFO")
		if err != nil {
			continue
		}
		iclnsig, err := row.Info().Get("CLNSIG")
		if err != nil {
			continue
		}
		mc := strings.Join(pkg.Interface2Array[string](imc), "|")
		geneinfos := strings.Split(strings.Join(pkg.Interface2Array[string](igeneinfo), "|"), "|")
		clnsig := strings.ToLower(strings.Join(pkg.Interface2Array[string](iclnsig), "|"))
		if strings.Contains(mc, "nonsense") || strings.Contains(mc, "frameshift") {
			is_pathogenic := strings.Contains(clnsig, "pathogenic") && !strings.Contains(clnsig, "conflict")
			for _, geneinfo := range geneinfos {
				geneid := strings.Split(geneinfo, ":")[1]
				_, ok := counts[geneid]
				if !ok {
					counts[geneid] = []int{0, 0}

				}
				counts[geneid][0] += 1
				if is_pathogenic {
					counts[geneid][1] += 1
				}
			}
		}
	}
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprint(writer, "GeneID\tTotal\tPathogenic\n")
	for geneid, count := range counts {
		fmt.Fprintf(writer, "%s\t%d\t%d\n", geneid, count[0], count[1])
	}
	return err
}

func NewPreClinvarGeneCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "gene",
		Short: "Prepare Clinvar Base on NCBI Clinvar VCF",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreClinvarGeneParam
			param.Input, _ = cmd.Flags().GetString("input")
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
	cmd.Flags().StringP("input", "i", "", "Input Clinvar VCF File")
	cmd.Flags().StringP("output", "o", "", "Output File")
	return cmd
}
