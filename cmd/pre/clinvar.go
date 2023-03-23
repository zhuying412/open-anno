package pre

import (
	"io"
	"log"
	"open-anno/pkg"
	"os"
	"path"

	"github.com/brentp/vcfgo"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type PreClinvarParam struct {
	Input  string `validate:"required,pathexists"`
	Output string `validate:"required"`
}

func (this PreClinvarParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this PreClinvarParam) HeaderInfoIDs() []string {
	return []string{"CLNREVSTAT", "CLNSIG", "CLNSIGCONF", "CLNDN"}
}

func (this PreClinvarParam) NewVcfWriter(writer io.WriteCloser, vcf string, infoKeys []string) (*vcfgo.Writer, error) {
	reader, err := pkg.NewIOReader(vcf)
	if err != nil {
		return &vcfgo.Writer{}, err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return &vcfgo.Writer{}, err
	}
	defer vcfReader.Close()
	vcfHeader := *vcfReader.Header
	vcfHeaderInfos := make(map[string]*vcfgo.Info)
	for key, info := range vcfHeader.Infos {
		if pkg.FindArr(infoKeys, key) != -1 {
			vcfHeaderInfos[info.Id] = &vcfgo.Info{
				Id:          info.Id,
				Description: info.Description,
				Type:        info.Type,
				Number:      info.Number,
			}
		}
	}
	vcfHeader.Infos = vcfHeaderInfos
	return vcfgo.NewWriter(writer, &vcfHeader)
}

func (this PreClinvarParam) Run() error {
	infoKeys := this.HeaderInfoIDs()
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
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	vcfWriter, err := this.NewVcfWriter(writer, this.Input, infoKeys)
	if err != nil {
		return err
	}
	for row := vcfReader.Read(); row != nil; row = vcfReader.Read() {
		row.Chromosome = "chr" + row.Chromosome
		for _, key := range row.Info().Keys() {
			if pkg.FindArr(infoKeys, key) == -1 {
				row.Info().Delete(key)
			}
		}
		vcfWriter.WriteVariant(row)
	}
	return err
}

func NewPreClinvarCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "cln",
		Short: "Prepare Clinvar Base on NCBI Clinvar VCF",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreClinvarParam
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
