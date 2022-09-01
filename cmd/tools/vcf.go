package tools

import (
	"errors"
	"log"
	"open-anno/pkg/io"

	"github.com/brentp/faidx"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type IAny2AnyParam interface {
	Valid() error
	Run() error
}

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
	TriosSamples []string `validate:"required"`
}

func (this TriosVCF2AnnoInputParam) Run() error {
	variants, err := io.ReadTriosVCF(this.Input, this.TriosSamples[0], this.TriosSamples[1], this.TriosSamples[2])
	if err != nil {
		return err
	}
	return io.WriteVariants(variants, this.Ouput)
}

type McVCF2AnnoInputParam struct {
	Any2AnyParam
	McSamples []string `validate:"required"`
}

func (this McVCF2AnnoInputParam) Run() error {
	variants, err := io.ReadMcVCF(this.Input, this.McSamples[0], this.McSamples[1])
	if err != nil {
		return err
	}
	return io.WriteVariants(variants, this.Ouput)
}

func NewVCF2AnnoInputCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "vcf2ai",
		Short: "Convert SNV VCF to AnnoInput",
		Run: func(cmd *cobra.Command, args []string) {
			triosSamples, _ := cmd.Flags().GetStringArray("trios")
			mcSamples, _ := cmd.Flags().GetStringArray("mother_child")
			var param IAny2AnyParam
			var baseParam Any2AnyParam
			baseParam.Input, _ = cmd.Flags().GetString("vcf")
			baseParam.Ouput, _ = cmd.Flags().GetString("avinput")
			if len(triosSamples) > 0 {
				param = TriosVCF2AnnoInputParam{Any2AnyParam: baseParam, TriosSamples: triosSamples}
			} else if len(mcSamples) > 0 {
				param = McVCF2AnnoInputParam{Any2AnyParam: baseParam, McSamples: mcSamples}
			} else {
				param = TriosVCF2AnnoInputParam{Any2AnyParam: baseParam}
			}
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
	cmd.Flags().StringArrayP("trios", "t", []string{}, "Trios Samples, the order must be Proband, Mother, Father")
	cmd.Flags().StringArrayP("mother_child", "m", []string{}, "Mother & Child Samples, the order must be Proband, Mother")
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
		Use:   "ai2vcf",
		Short: "Convert SNV AnnoInput to VCF",
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
