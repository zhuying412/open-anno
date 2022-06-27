package pre

import (
	"bufio"
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"strconv"
	"strings"

	"github.com/spf13/cobra"
)

func RunIndexDatabase(infile string, binSize int) {
	log.Println("Init parameters ...")
	idxfile := infile + ".idx"
	reader, err := io.NewIoReader(infile)
	if err != nil {
		log.Fatal(err)
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
			log.Fatal(err)
		}
		curbin := pkg.CurBin(chrom, start, binSize)
		if _, ok := idxMap[curbin]; ok {
			idxMap[curbin].End = offset + length
		} else {
			idxs = append(idxs, curbin)
			idxMap[curbin] = &io.DBVarIdx{Bin: curbin, Start: offset, End: offset + length}
		}
		offset += length
	}
	writer, err := io.NewIoWriter(idxfile)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintf(writer, "#Bin\t%d\n", binSize)
	for _, key := range idxs {
		idx := idxMap[key]
		fmt.Fprintf(writer, "%s\t%d\t%d\n", idx.Bin, idx.Start, idx.End)
	}
}

func NewIndexDatabaseCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "idx",
		Short: "Index FilterBased or RegionBased database",
		Run: func(cmd *cobra.Command, args []string) {
			infile, _ := cmd.Flags().GetString("infile")
			binSize, _ := cmd.Flags().GetInt("binsize")
			if infile == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				RunIndexDatabase(infile, binSize)
			}
		},
	}
	cmd.Flags().StringP("infile", "i", "", "Input Fileter-Based File")
	cmd.Flags().IntP("binsize", "b", 1000, "Bin Size")
	return cmd
}
