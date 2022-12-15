package pre

import (
	"fmt"
	"log"
	"open-anno/anno"
	"open-anno/pkg"
	"os"
	"path"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type PreGenePredParam struct {
	Input     string `validate:"required,pathexists"`
	RefDict   string `validate:"required,pathexists"`
	Output    string `validate:"required"`
	IndexStep int    `validate:"required"`
}

func (this PreGenePredParam) OutIdx() string {
	return this.Output + ".idx"
}

func (this PreGenePredParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this PreGenePredParam) CreateGenePred(genome map[string]int) error {
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
		row := strings.Split(scanner.Text(), "\t")
		if _, ok := genome[row[2]]; ok {
			fmt.Fprintln(writer, scanner.Text())
		}
	}
	return nil
}

func (this PreGenePredParam) Run() error {
	// 读取RefDict
	log.Printf("Read RefDict: %s ...", this.RefDict)
	genome, err := anno.ReadGenomeDict(this.RefDict)
	if err != err {
		return err
	}
	// 创建 GenePred
	if this.Input != this.Output {
		log.Printf("Create GenePred: %s ...", this.Output)
		this.CreateGenePred(genome)
	}
	// 构建GenePred索引
	log.Printf("Read GenePred: %s ...", this.Input)
	gpes, err := pkg.ReadGenePred(this.Input)
	if err != err {
		return err
	}
	log.Printf("Create TransIndex: %s ...", this.OutIdx())
	return pkg.CreateTransIndexes(gpes, genome, this.IndexStep, this.OutIdx())
}

func NewPreGenePredCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "gpe",
		Short: "Prepare GenePred Files",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreGenePredParam
			param.RefDict, _ = cmd.Flags().GetString("refdict")
			param.Input, _ = cmd.Flags().GetString("input")
			param.Output, _ = cmd.Flags().GetString("output")
			param.IndexStep, _ = cmd.Flags().GetInt("step")
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
	cmd.Flags().StringP("refdict", "r", "", "Input Reference Dict File")
	cmd.Flags().StringP("input", "i", "", "Input GenePred File")
	cmd.Flags().StringP("output", "o", "", "Output GenePred File")
	cmd.Flags().IntP("step", "l", 300000, "Transcript Index Step Length")
	return cmd
}
