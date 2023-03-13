package pre

import (
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"open-anno/pkg"
	"os"
	"path"
	"strings"

	"github.com/brentp/vcfgo"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type PreGnomadParam struct {
	Input  string `validate:"required,pathexists"`
	Output string `validate:"required"`
}

func (this PreGnomadParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this PreGnomadParam) Inputs() ([]string, error) {
	vcfs := make([]string, 0)
	fileinfos, err := ioutil.ReadDir(this.Input)
	if err != nil {
		return vcfs, err
	}
	for _, file := range fileinfos {
		if file.IsDir() != true && strings.HasSuffix(file.Name(), "vcf.bgz") {
			vcfs = append(vcfs, path.Join(this.Input, file.Name()))
		}
	}
	return vcfs, nil
}

func (this PreGnomadParam) HeaderInfoIDs() []string {
	populations := []string{"", "oth", "ami", "sas", "fin", "eas", "amr", "afr", "asj", "nfe"}
	prefixes := []string{"AC", "AN", "AF", "nhomalt"}
	keys := make([]string, len(populations)*len(prefixes))
	for i, population := range populations {
		for j, prefix := range prefixes {
			key := prefix
			if population != "" {
				key += "_" + population
			}
			keys[len(prefixes)*i+j] = key
		}
	}
	return keys
}

func (this PreGnomadParam) NewVcfWriter(writer io.WriteCloser, vcf string, infoKeys []string, dbname string) (*vcfgo.Writer, error) {
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
			id := dbname + "_" + key
			vcfHeaderInfos[id] = &vcfgo.Info{
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

func (this PreGnomadParam) Run() error {
	dbname := "gnomAD"
	vcfs, err := this.Inputs()
	if err != nil {
		return err
	}
	infoKeys := this.HeaderInfoIDs()
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	vcfWriter, err := this.NewVcfWriter(writer, vcfs[0], infoKeys, dbname)
	if err != nil {
		return err
	}
	for _, vcf := range vcfs {
		fmt.Println(vcf)
		reader, err := pkg.NewIOReader(vcf)
		if err != nil {
			return err
		}
		// defer reader.Close()
		vcfReader, err := vcfgo.NewReader(reader, false)
		if err != nil {
			return err
		}
		for row := vcfReader.Read(); row != nil; row = vcfReader.Read() {
			for _, key := range row.Info().Keys() {
				if pkg.FindArr(infoKeys, key) == -1 {
					row.Info().Delete(key)
				}
			}
			vcfWriter.WriteVariant(row)
		}
	}
	return nil
}

func NewPreGnomadCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "gnomad",
		Short: "Prepare GnomAD Frequency Base on GnomAD VCF",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreGnomadParam
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
	cmd.Flags().StringP("input", "i", "", "Input gnomAD VCF Directory")
	cmd.Flags().StringP("output", "o", "", "Output File")
	return cmd
}
