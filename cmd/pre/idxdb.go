package pre

import (
	"fmt"
	"log"
	"open-anno/anno/db"
	"open-anno/pkg"
	"strconv"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type IdxDBParam struct {
	Input   string `validate:"required,pathexists"`
	BinSize int    `validate:"required"`
}

func (this IdxDBParam) OutIdx() string {
	return this.Input + ".idx"
}

func (this IdxDBParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	return validate.Struct(this)
}

func (this IdxDBParam) Run() error {
	log.Println("Start indexing database file")
	reader, err := pkg.NewIOReader(this.Input)
	if err != nil {
		return err
	}
	defer reader.Close()
	var offset int64
	idxMap := make(map[string]*db.FilterVarIdx)
	indexes := make([]string, 0)
	scanner := pkg.NewIOScanner(reader)
	var lastChrom string
	for scanner.Scan() {
		line := scanner.Text()
		length := int64(len(line) + 1)
		if strings.HasPrefix(line, "#") {
			offset += length
			continue
		}
		row := strings.Split(line, "\t")
		chrom := row[0]
		if chrom != lastChrom {
			fmt.Printf("Reading %s\n", chrom)
			lastChrom = chrom
		}
		start, err := strconv.Atoi(row[1])
		if err != nil {
			return err
		}
		curbin := db.CurBin(chrom, start, this.BinSize)
		if _, ok := idxMap[curbin]; !ok {
			indexes = append(indexes, curbin)
			idxMap[curbin] = &db.FilterVarIdx{Bin: curbin, Start: offset}
		}
		idxMap[curbin].End = offset + length
		offset += length
	}
	writer, err := pkg.NewIOWriter(this.OutIdx())
	if err != nil {
		return err
	}
	fmt.Fprintf(writer, "#Bin\t%d\n", this.BinSize)
	for _, key := range indexes {
		idx := idxMap[key]
		fmt.Fprintf(writer, "%s\t%d\t%d\n", idx.Bin, idx.Start, idx.End)
	}
	return err
}

func NewIdxDBCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "idx",
		Short: "Index FilterBased or RegionBased database",
		Run: func(cmd *cobra.Command, args []string) {
			var param IdxDBParam
			param.Input, _ = cmd.Flags().GetString("infile")
			param.BinSize, _ = cmd.Flags().GetInt("binsize")
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
	cmd.Flags().StringP("infile", "i", "", "Input FilterBased File")
	cmd.Flags().IntP("binsize", "b", 1000, "Bin Size")
	return cmd
}
