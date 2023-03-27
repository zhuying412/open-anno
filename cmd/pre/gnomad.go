package pre

import (
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
	return os.MkdirAll(this.Output, 0666)
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
	prefixes := []string{"AC", "AN", "AF", "nhomalt"}
	populations := []string{"", "oth", "ami", "sas", "fin", "eas", "amr", "afr", "asj", "nfe"}
	sexes := []string{"", "XY"}
	keys := make([]string, 0)
	for _, prefix := range prefixes {
		for _, population := range populations {
			for _, sex := range sexes {
				key := prefix
				if population != "" {
					key += "_" + population
				}
				if sex != "" {
					key += "_" + sex
				}
				keys = append(keys, key)
			}
		}
	}
	return keys
}

func (this PreGnomadParam) ProcessVCF(inVcf string, outVcf string, infoKeys []string, dbname string, errChan chan error) {
	reader, err := pkg.NewIOReader(inVcf)
	if err != nil {
		errChan <- err
		return
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		errChan <- err
		return
	}
	defer vcfReader.Close()
	vcfHeaderInfos := make(map[string]*vcfgo.Info)
	for key, info := range vcfReader.Header.Infos {
		if pkg.FindArr(infoKeys, key) != -1 {
			id := dbname + "_" + key
			vcfHeaderInfos[id] = &vcfgo.Info{
				Id:          id,
				Description: info.Description,
				Type:        info.Type,
				Number:      info.Number,
			}
		}
	}

	vcfHeader := &vcfgo.Header{
		SampleNames:   vcfReader.Header.SampleNames,
		SampleFormats: vcfReader.Header.SampleFormats,
		Filters:       vcfReader.Header.Filters,
		Extras:        vcfReader.Header.Extras,
		FileFormat:    vcfReader.Header.FileFormat,
		Contigs:       vcfReader.Header.Contigs,
		Samples:       vcfReader.Header.Samples,
		Pedigrees:     vcfReader.Header.Pedigrees,
		Infos:         vcfHeaderInfos,
	}
	writer, err := pkg.NewIOWriter(outVcf)
	if err != nil {
		errChan <- err
		return
	}
	vcfWriter, err := vcfgo.NewWriter(writer, vcfHeader)
	if err != nil {
		errChan <- err
		return
	}
	for row := vcfReader.Read(); row != nil; row = vcfReader.Read() {
		variant := &vcfgo.Variant{
			Chromosome: row.Chromosome,
			Pos:        row.Pos,
			Reference:  row.Reference,
			Alternate:  row.Alternate,
			Info_:      &vcfgo.InfoByte{},
		}
		for _, key := range infoKeys {
			val, err := row.Info().Get(key)
			if err != nil {
				continue
			}
			variant.Info().Set(dbname+"_"+key, val)
		}
		vcfWriter.WriteVariant(variant)
	}
	errChan <- nil
	return
}

func (this PreGnomadParam) Run() error {
	dbname := "gnomAD"
	vcfs, err := this.Inputs()
	if err != nil {
		return err
	}
	infoKeys := this.HeaderInfoIDs()
	errChan := make(chan error, len(vcfs))
	for _, vcf := range vcfs {
		go this.ProcessVCF(vcf, path.Join(this.Output, strings.ReplaceAll(path.Base(vcf), ".bgz", "")), infoKeys, dbname, errChan)
	}
	for i := 0; i < len(vcfs); i++ {
		err = <-errChan
		if err != nil {
			return err
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
	cmd.Flags().StringP("output", "o", "", "Output Directory")
	return cmd
}
