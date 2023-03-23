package pre

// import (
// 	"fmt"
// 	"log"
// 	"open-anno/pkg"
// 	"os"
// 	"path"

// 	"github.com/brentp/vcfgo"
// 	"github.com/go-playground/validator/v10"
// 	"github.com/spf13/cobra"
// )

// type Pre1000GenomesParam struct {
// 	Input  string `validate:"required,pathexists"`
// 	Output string `validate:"required"`
// }

// func (this Pre1000GenomesParam) Valid() error {
// 	validate := validator.New()
// 	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
// 	err := validate.Struct(this)
// 	if err != nil {
// 		return err
// 	}
// 	outdir := path.Dir(this.Output)
// 	return os.MkdirAll(outdir, 0666)
// }

// func (this Pre1000GenomesParam) Run() error {
// 	reader, err := pkg.NewIOReader(this.Input)
// 	if err != nil {
// 		return err
// 	}
// 	defer reader.Close()
// 	vcfReader, err := vcfgo.NewReader(reader, false)
// 	if err != nil {
// 		return err
// 	}
// 	defer vcfReader.Close()
// 	writer, err := pkg.NewIOWriter(this.Output)
// 	if err != nil {
// 		return err
// 	}
// 	defer writer.Close()
// 	fmt.Fprint(writer, "#Chr\tStart\tEnd\tRef\tAlt\t1000Genomes_AC\t1000Genomes_AN\t1000Genomes_AF\t1000Genomes_AF_EAS\t1000Genomes_AF_EUR\t1000Genomes_AF_AFR\t1000Genomes_AF_AMR\t1000Genomes_AF_SAS\n")
// 	for {
// 		row := vcfReader.Read()
// 		if row == nil {
// 			break
// 		}
// 		for i, alt := range row.Alt() {
// 			chrom, start, end, ref, alt := pkg.VCFtoAV(row.Chrom(), int(row.Pos), row.Ref(), alt)
// 			ac, err := row.Info().Get("AC")
// 			if err != nil {
// 				return err
// 			}
// 			if v, ok := ac.([]int); ok {
// 				ac = v[i]
// 			}
// 			an, err := row.Info().Get("AN")
// 			if err != nil {
// 				return err
// 			}
// 			if v, ok := an.([]int); ok {
// 				an = v[i]
// 			}
// 			af, err := row.Info().Get("AF")
// 			if err != nil {
// 				return err
// 			}
// 			if v, ok := af.([]float32); ok {
// 				af = v[i]
// 			}
// 			afEAS, err := row.Info().Get("EAS_AF")
// 			if err != nil {
// 				return err
// 			}
// 			if v, ok := afEAS.([]float32); ok {
// 				afEAS = v[i]
// 			}
// 			afEUR, err := row.Info().Get("EUR_AF")
// 			if err != nil {
// 				return err
// 			}
// 			if v, ok := afEUR.([]float32); ok {
// 				afEUR = v[i]
// 			}
// 			afAFR, err := row.Info().Get("AFR_AF")
// 			if err != nil {
// 				return err
// 			}
// 			if v, ok := afAFR.([]float32); ok {
// 				afAFR = v[i]
// 			}
// 			afAMR, err := row.Info().Get("AMR_AF")
// 			if err != nil {
// 				return err
// 			}
// 			if v, ok := afAMR.([]float32); ok {
// 				afAMR = v[i]
// 			}
// 			afSAS, err := row.Info().Get("SAS_AF")
// 			if err != nil {
// 				return err
// 			}
// 			if v, ok := afSAS.([]float32); ok {
// 				afSAS = v[i]
// 			}
// 			fmt.Fprintf(writer,
// 				"chr%s\t%d\t%d\t%s\t%s\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
// 				chrom, start, end, ref, alt, ac, an, af,
// 				afEAS, afEUR, afAFR, afAMR, afSAS,
// 			)
// 		}
// 	}
// 	return err
// }

// func NewPre1000GenomesCmd() *cobra.Command {
// 	cmd := &cobra.Command{
// 		Use:   "1kg",
// 		Short: "Prepare 1000 Genomes Frequency Base on 1000 Genomes VCF",
// 		Run: func(cmd *cobra.Command, args []string) {
// 			var param Pre1000GenomesParam
// 			param.Input, _ = cmd.Flags().GetString("input")
// 			param.Output, _ = cmd.Flags().GetString("output")
// 			err := param.Valid()
// 			if err != nil {
// 				cmd.Help()
// 				log.Fatal(err)
// 			}
// 			err = param.Run()
// 			if err != nil {
// 				log.Fatal(err)
// 			}
// 		},
// 	}
// 	cmd.Flags().StringP("input", "i", "", "Input 1000Genomes VCF File")
// 	cmd.Flags().StringP("output", "o", "", "Output File")
// 	return cmd
// }
