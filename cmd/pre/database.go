package pre

import (
	"bufio"
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/gene"
	"os"
	"path"
	"strings"

	"github.com/spf13/cobra"
)

func RunPreDatabase(infile string, builder string, outdir string) {
	log.Println("Init parameters ...")
	gene.SetGenome(builder)
	if _, err := os.Stat(outdir); os.IsNotExist(err) {
		err := os.MkdirAll(outdir, os.ModePerm)
		if err != nil {
			log.Fatal(err)
		}
	}
	log.Printf("Write database: %s ...", outdir)
	fi, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer fi.Close()
	scanner := bufio.NewScanner(fi)
	scanner.Scan()
	header := scanner.Text()
	if !strings.HasPrefix(header, "#Chr") {
		log.Fatalf("error database file, header not found: %s", infile)
	}
	writers := make(map[string]*os.File)
	for scanner.Scan() {
		line := scanner.Text()
		fileds := strings.Split(line, "\t")
		chrom := pkg.FormatChrom(fileds[0])
		if _, ok := gene.GENOME[chrom]; ok {
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

func NewPreDatabaseCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "db",
		Short: "Prepare required FilterBased or RegionBased database",
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
				RunPreDatabase(infile, strings.ToLower(builder), path.Join(dbpath, builder, name))
			}
		},
	}
	cmd.Flags().StringP("infile", "i", "", "Input Fileter-Based File")
	cmd.Flags().StringP("dbpath", "d", "", "Database Directory")
	cmd.Flags().StringP("name", "n", "", "Database Name")
	cmd.Flags().StringP("builder", "b", "hg19", "Database Path")
	return cmd
}
