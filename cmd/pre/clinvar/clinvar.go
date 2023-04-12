package clinvar

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"os"
	"path"
	"regexp"
	"strings"

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
	return []string{"ALLELEID", "CLNREVSTAT", "CLNSIG", "CLNSIGCONF", "CLNDN"}
}

func (this PreClinvarParam) Run() error {
	infoKeys := this.HeaderInfoIDs()
	reader, err := pkg.NewIOReader(this.Input)
	if err != nil {
		return err
	}
	defer reader.Close()
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		text := scanner.Text()
		if strings.HasPrefix(text, "#") {
			if strings.HasPrefix(text, "##INFO=") {
				for _, infoKey := range infoKeys {
					if strings.HasPrefix(text, "##INFO=<ID="+infoKey) {
						fmt.Fprintln(writer, strings.Replace(text, infoKey, "ClinVar_"+infoKey, 1))
						break
					}
				}
			} else {
				fmt.Fprintln(writer, text)
			}
		} else {
			row := strings.Split(text, "\t")
			infos := make([]string, 0)
			for _, item := range strings.Split(row[7], ";") {
				info := strings.Split(item, "=")
				if pkg.FindArr(infoKeys, info[0]) != -1 {
					info[1] = "ClinVar_" + regexp.MustCompile(`[^\w\(\)]+`).ReplaceAllString(info[1], "|")
					infos = append(infos, fmt.Sprintf("ClinVar_%s=%s", info[0], info[1]))
				}
			}
			row[7] = strings.Join(infos, ";")
			fmt.Fprintln(writer, strings.Join(row, "\t"))
		}
	}
	return err
}

func NewPreClinvarCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "db",
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
