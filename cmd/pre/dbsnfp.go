package pre

import (
	"log"
	"open-anno/pkg"
	"os"
	"path"
	"strconv"
	"strings"

	"github.com/brentp/vcfgo"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type PreDbnsfpParam struct {
	Input  string `validate:"required,pathexists"`
	Output string `validate:"required"`
}

func (this PreDbnsfpParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this PreDbnsfpParam) HeaderInfos(fieldNames []string) map[string]*vcfgo.Info {
	headerInfos := make(map[string]*vcfgo.Info)
	for _, name := range fieldNames[5:] {
		info := vcfgo.Info{
			Id:          name,
			Number:      "A",
			Type:        "Float",
			Description: name,
		}
		if strings.HasSuffix(name, "pred") || strings.HasPrefix(name, "GTEx_V8") || name == "Interpro_domain" {
			info.Number = "."
			info.Type = "String"
		}
		headerInfos[name] = &info
	}
	return headerInfos
}

func (this PreDbnsfpParam) Run() error {
	reader, err := pkg.NewIOReader(this.Input)
	if err != nil {
		return err
	}
	defer reader.Close()
	scanner := pkg.NewCSVScanner(reader)
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	vcfWriter, err := vcfgo.NewWriter(writer, &vcfgo.Header{Infos: this.HeaderInfos(scanner.FieldNames), FileFormat: "4.1"})
	if err != nil {
		return err
	}
	for scanner.Scan() {
		row := scanner.Row()
		pos, err := strconv.Atoi(row["Start"])
		if err != nil {
			return err
		}
		variant := &vcfgo.Variant{
			Chromosome: row["#Chr"],
			Pos:        uint64(pos),
			Id_:        ".",
			Reference:  row["Ref"],
			Alternate:  []string{row["Alt"]},
			Info_:      &vcfgo.InfoByte{},
		}
		for _, name := range scanner.FieldNames[5:] {
			val := row[name]
			if val[0] != '.' {
				variant.Info().Set(name, strings.ReplaceAll(val, ";", "|"))
			}
		}
		vcfWriter.WriteVariant(variant)
	}
	return nil
}

func NewPreDbnsfpCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "dbnsfp",
		Short: "Prepare dbNSFP VCF",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreDbnsfpParam
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
	cmd.Flags().StringP("input", "i", "", "Input ANNOVAR dbNSFP Database File")
	cmd.Flags().StringP("output", "o", "", "Output VCF File")
	return cmd
}
