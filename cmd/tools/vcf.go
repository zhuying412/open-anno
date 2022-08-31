package tools

import (
	"errors"
	"log"
	"open-anno/pkg/io"

	"github.com/brentp/faidx"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type Any2AnyParam struct {
	Input string `validate:"required"`
	Ouput string `validate:"required"`
}

func (this Any2AnyParam) Valid() error {
	validate := validator.New()
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return nil
}

func Run() error {
	return errors.New("Not Impl")
}

type VCF2AnnoInputParam struct {
	Any2AnyParam
}

func (this VCF2AnnoInputParam) Run() error {
	variants, err := io.ReadVCF(this.Input)
	if err != nil {
		return err
	}
	return io.WriteVariants(variants, this.Ouput)
}

type TriosVCF2AnnoInputParam struct {
	Any2AnyParam
	Proband string `validate:"required"`
	Mother  string `validate:"required"`
	Father  string `validate:"required"`
}

func (this TriosVCF2AnnoInputParam) Run() error {
	variants, err := io.ReadTriosVCF(this.Input, this.Proband, this.Mother, this.Father)
	if err != nil {
		return err
	}
	return io.WriteVariants(variants, this.Ouput)
}

func NewVCF2AnnoInputCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "vcf2av",
		Short: "Convert SNV VCF to AVINPUT",
		Run: func(cmd *cobra.Command, args []string) {
			is_trios, _ := cmd.Flags().GetBool("trios")
			if is_trios {
				var param TriosVCF2AnnoInputParam
				param.Input, _ = cmd.Flags().GetString("vcf")
				param.Ouput, _ = cmd.Flags().GetString("avinput")
				param.Proband, _ = cmd.Flags().GetString("mother")
				param.Father, _ = cmd.Flags().GetString("father")
				param.Mother, _ = cmd.Flags().GetString("mother")
				err := param.Valid()
				if err != nil {
					cmd.Help()
					log.Fatal(err)
				}
				err = param.Run()
				if err != nil {
					log.Fatal(err)
				}
			} else {
				var param VCF2AnnoInputParam
				param.Input, _ = cmd.Flags().GetString("vcf")
				param.Ouput, _ = cmd.Flags().GetString("avinput")
				err := param.Valid()
				if err != nil {
					cmd.Help()
					log.Fatal(err)
				}
				err = param.Run()
				if err != nil {
					log.Fatal(err)
				}
			}

		},
	}
	cmd.Flags().StringP("vcf", "i", "", "VCF File")
	cmd.Flags().StringP("avinput", "o", "", "AVINPUT File")
	cmd.Flags().BoolP("trios", "t", false, "Is Trios VCF or not")
	cmd.Flags().StringP("proband", "p", "", "Proband Sample Name, -t required")
	cmd.Flags().StringP("father", "f", "", "Father Sample Name, -t required")
	cmd.Flags().StringP("mother", "m", "", "Mother Sample Name, -t required")
	return cmd
}

type Av2VcfParam struct {
	Any2AnyParam
	Genome string `validate:"required"`
}

func (this Av2VcfParam) Run() error {
	fai, err := faidx.New(this.Genome)
	if err != nil {
		return err
	}
	variants, err := io.ReadVariants(this.Input, this.Ouput)
	if err != nil {
		return err
	}
	return io.WriteVariantsToVCF(variants, fai, this.Ouput)
}

func NewAnnoInput2VCFCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "av2vcf",
		Short: "Convert SNV AVINPUT to VCF",
		Run: func(cmd *cobra.Command, args []string) {
			var param Av2VcfParam
			param.Input, _ = cmd.Flags().GetString("avinput")
			param.Ouput, _ = cmd.Flags().GetString("vcf")
			param.Genome, _ = cmd.Flags().GetString("genome")
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
	cmd.Flags().StringP("vcf", "i", "", "VCF File")
	cmd.Flags().StringP("avinput", "o", "", "AVINPUT File")
	cmd.Flags().StringP("genome", "g", "", "Genome Fasta File")
	return cmd
}
