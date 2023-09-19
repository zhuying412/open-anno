package tools

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"os"
	"path"
	"sort"
	"strconv"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type Transcript struct {
	Chrom  string
	Gene   string
	Name   string
	MANE   int
	HGMD   int
	Length int
}

func (this Transcript) Source() string {
	if this.MANE == 2 {
		return "MANE_Plus_Clinical"
	} else if this.MANE == 1 {
		return "MANE_Select"
	} else {
		if this.HGMD == 1 {
			return "HGMD"
		}
	}
	return "Max_Length"
}

type Transcripts []Transcript

func (this Transcripts) Len() int      { return len(this) }
func (this Transcripts) Swap(i, j int) { this[i], this[j] = this[j], this[i] }
func (this Transcripts) Less(i, j int) bool {
	if this[i].MANE == this[j].MANE {
		if this[i].HGMD == this[j].HGMD {
			return this[i].Length > this[j].Length
		}
		return this[i].HGMD > this[j].HGMD
	}
	return this[i].MANE > this[j].MANE
}

type RepTransParam struct {
	RefseqCurated string `validate:"required,pathexists"`
	MANE          string `validate:"required,pathexists"`
	RefseqHGMD    string `validate:"required,pathexists"`
	Output        string `validate:"required"`
}

func (this RepTransParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this RepTransParam) ReadMANE(maneFile string) (map[string]int, error) {
	data := make(map[string]int)
	reader, err := pkg.NewIOReader(maneFile)
	if err != nil {
		return data, err
	}
	defer reader.Close()
	scanner := pkg.NewCSVScanner(reader)
	for scanner.Scan() {
		row := scanner.Row()
		chrom := row["GRCh38_chr"]
		if strings.HasPrefix(chrom, "NC_") {
			chrom = strings.Split(chrom, ".")[0][7:9]
			if chrom[0] == '0' {
				chrom = chrom[1:]
			}
			if chrom == "23" {
				chrom = "X"
			}
			if chrom == "24" {
				chrom = "Y"
			}
			chrom = "chr" + chrom
			gene := row["symbol"]
			trans := strings.Split(row["RefSeq_nuc"], ".")[0]
			key := fmt.Sprintf("%s\t%s\t%s", chrom, gene, trans)
			val := 1
			if row["MANE_status"] == "MANE Plus Clinical" {
				val = 2
			}
			data[key] = pkg.Max(data[key], val)
		}
	}
	return data, nil
}

func (this RepTransParam) ReadHGMD(hgmdGpeFile string) (map[string]int, error) {
	data := make(map[string]int)
	reader, err := pkg.NewIOReader(hgmdGpeFile)
	if err != nil {
		return data, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		chrom := row[2]
		if len(chrom) <= 5 {
			gene := row[12]
			trans := strings.Split(row[1], ".")[0]
			key := fmt.Sprintf("%s\t%s\t%s", chrom, gene, trans)
			data[key] = 1
		}
	}
	return data, nil
}

func (this RepTransParam) ReadCuratedRefSeq(gpeFile string, maneData, hgmdData map[string]int) (map[string]Transcripts, error) {
	data := make(map[string]Transcripts)
	reader, err := pkg.NewIOReader(gpeFile)
	if err != nil {
		return data, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		chrom := row[2]
		if len(chrom) <= 5 {
			gene := row[12]
			name := strings.Split(row[1], ".")[0]
			start, err := strconv.Atoi(row[4])
			if err != nil {
				return data, err
			}
			end, err := strconv.Atoi(row[5])
			if err != nil {
				return data, err
			}
			key := fmt.Sprintf("%s\t%s\t%s", chrom, gene, name)
			trans := Transcript{
				Chrom:  chrom,
				Gene:   gene,
				Name:   row[1],
				Length: end - start,
				MANE:   maneData[key],
				HGMD:   hgmdData[key],
			}
			key2 := fmt.Sprintf("%s\t%s", chrom, gene)
			if _, ok := data[key2]; !ok {
				data[key2] = Transcripts{trans}
			} else {
				data[key2] = append(data[key2], trans)
			}
		}
	}
	return data, nil
}

func (this RepTransParam) Run() error {
	maneData, err := this.ReadMANE(this.MANE)
	if err != nil {
		return err
	}
	hgmdData, err := this.ReadHGMD(this.RefseqHGMD)
	if err != nil {
		return err
	}
	gpeData, err := this.ReadCuratedRefSeq(this.RefseqCurated, maneData, hgmdData)
	if err != nil {
		return err
	}
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	for _, transcripts := range gpeData {
		sort.Sort(transcripts)
		trans := transcripts[0]
		fmt.Fprintf(writer, "%s\t%s\t%s\t%s\n", trans.Chrom, trans.Gene, trans.Name, trans.Source())
	}
	return nil
}

func NewRepTransCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "rt",
		Short: "Get Representative Transcript",
		Run: func(cmd *cobra.Command, args []string) {
			var param RepTransParam
			param.RefseqCurated, _ = cmd.Flags().GetString("curated")
			param.RefseqHGMD, _ = cmd.Flags().GetString("hgmd")
			param.MANE, _ = cmd.Flags().GetString("mane")
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
	cmd.Flags().StringP("curated", "c", "", "Input Refseq Curated File from UCSC, name: ncbiRefSeqCurated.txt.gz")
	cmd.Flags().StringP("mane", "m", "", "Input Refseq Select File from UCSC, name: ncbiRefSeqSelect.txt.gz")
	cmd.Flags().StringP("hgmd", "g", "", "Input Refseq HGMD File  from UCSC, name: ncbiRefSeqHgmd.txt.gz")
	cmd.Flags().StringP("output", "o", "-", "Output Representative Transcript File")
	return cmd
}
