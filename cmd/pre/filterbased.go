package pre

import (
	"bufio"
	"fmt"
	"log"
	"open-anno/pkg/gene"
	"os"
	"path"
	"strings"

	"github.com/spf13/cobra"
)

func PreFilterBased(infile string, builder string, outdir string) {
	log.Println("Init parameters ...")
	genome := gene.NewGenome(builder)
	if _, err := os.Stat(outdir); os.IsNotExist(err) {
		err := os.MkdirAll(outdir, os.ModePerm)
		if err != nil {
			log.Fatal(err)
		}
	}
	log.Printf("Write filterbased: %s ...", outdir)
	fi, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer fi.Close()
	scanner := bufio.NewScanner(fi)
	writers := make(map[string]*os.File)
	var header string
	for scanner.Scan() {
		line := scanner.Text()
		fileds := strings.Split(line, "\t")
		chrom := fileds[0]
		if chrom == "Chr" {
			header = line
			continue
		}
		if _, ok := genome[chrom]; ok {
			if _, ok := writers[chrom]; !ok {
				outfile := path.Join(outdir, fmt.Sprintf("chr%s.txt", chrom))
				writers[chrom], err = os.Create(outfile)
				if err != nil {
					log.Fatal(err)
				}
				writers[chrom].WriteString(header + "\n")
			}
			writers[chrom].WriteString(line + "\n")
		}

	}
	for _, writer := range writers {
		err := writer.Close()
		if err != nil {
			log.Panic()
		}
	}
}

func NewPreFilterBasedCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "FitlerBased",
		Short: "Prepare required FilterBased database",
		Run: func(cmd *cobra.Command, args []string) {
			infile, _ := cmd.Flags().GetString("infile")
			dbpath, _ := cmd.Flags().GetString("dbpath")
			builder, _ := cmd.Flags().GetString("builder")
			name, _ := cmd.Flags().GetString("name")
			if infile == "" || dbpath == "" || builder == "" || name == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				PreFilterBased(infile, builder, path.Join(dbpath, builder, name))
			}
		},
	}
	cmd.Flags().StringP("infile", "i", "", "Input Fileter-Based File")
	cmd.Flags().StringP("dbpath", "d", "", "Database Directory")
	cmd.Flags().StringP("name", "n", "", "Database Name")
	cmd.Flags().StringP("builder", "b", "hg19", "Database Path")
	return cmd
}
