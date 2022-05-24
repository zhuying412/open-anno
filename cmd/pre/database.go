package pre

import (
	"bufio"
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"open-anno/pkg/seq"
	"os"
	"path"
	"strings"

	"github.com/spf13/cobra"
)

func RunPreDatabase(infile string, builder string, outdir string) {
	log.Println("Init parameters ...")
	seq.SetGenome(builder)
	if _, err := os.Stat(outdir); os.IsNotExist(err) {
		err := os.MkdirAll(outdir, os.ModePerm)
		if err != nil {
			log.Fatal(err)
		}
	}
	log.Printf("Write database: %s ...", outdir)
	reader, err := io.NewIoReader(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer reader.Close()
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	scanner.Scan()
	header := scanner.Text()
	if !strings.HasPrefix(header, "#Chr") {
		log.Fatalf("error database file, header not found: %s", infile)
	}
	writers := make(map[string]io.WriteCloser)
	for scanner.Scan() {
		line := scanner.Text()
		fileds := strings.Split(line, "\t")
		chrom := pkg.FormatChrom(fileds[0])
		if _, ok := seq.GENOME[chrom]; ok {
			if _, ok := writers[chrom]; !ok {
				outfile := path.Join(outdir, fmt.Sprintf("chr%s.txt.gz", chrom))
				writers[chrom], err = io.NewIoWriter(outfile)
				fmt.Fprintf(writers[chrom], "%s\n", header)
			}
			fmt.Fprintf(writers[chrom], "%s\n", line)
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
