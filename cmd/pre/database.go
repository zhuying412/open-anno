package pre

import (
	"bufio"
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"strconv"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type IdxDBParam struct {
	Input   string `validate:"required,pathexists"`
	BinSize int    `validate:"required"`
}

func (this IdxDBParam) OutIndex() string {
	return this.Input + ".idx"
}

func (this IdxDBParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return nil
}

func (this IdxDBParam) Run() error {
	log.Println("Init parameters ...")
	reader, err := io.NewIoReader(this.Input)
	if err != nil {
		return err
	}
	defer reader.Close()
	var offset int64
	idxMap := make(map[string]*io.DBVarIdx)
	idxs := make([]string, 0)
	scanner := bufio.NewScanner(reader)
	for scanner.Scan() {
		line := scanner.Text()
		length := int64(len(line) + 1)
		if strings.HasPrefix(line, "#") {
			offset += length
			continue
		}
		field := strings.Split(line, "\t")
		chrom := field[0]
		start, err := strconv.Atoi(field[1])
		if err != nil {
			return err
		}
		curbin := pkg.CurBin(chrom, start, this.BinSize)
		if _, ok := idxMap[curbin]; ok {
			idxMap[curbin].End = offset + length
		} else {
			idxs = append(idxs, curbin)
			idxMap[curbin] = &io.DBVarIdx{Bin: curbin, Start: offset, End: offset + length}
		}
		offset += length
	}
	writer, err := io.NewIoWriter(this.OutIndex())
	if err != nil {
		return err
	}
	fmt.Fprintf(writer, "#Bin\t%d\n", this.BinSize)
	for _, key := range idxs {
		idx := idxMap[key]
		fmt.Fprintf(writer, "%s\t%d\t%d\n", idx.Bin, idx.Start, idx.End)
	}
	return err
}

func NewIndexDatabaseCmd() *cobra.Command {
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
	cmd.Flags().StringP("infile", "i", "", "Input Fileter-Based File")
	cmd.Flags().IntP("binsize", "b", 1000, "Bin Size")
	return cmd
}
