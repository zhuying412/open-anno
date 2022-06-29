package cmd

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

func IndexDatabse(infile string, bin int) error {
	reader, err := io.NewIoReader(infile)
	if err != nil {
		return err
	}
	defer reader.Close()
	var offset int
	binMap := make(map[string][]int)
	bins := make([]string, 0)
	rd := bufio.NewReader(reader)
	for {
		line, err := rd.ReadString('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		length := len(line)
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
		curbin := pkg.CurBin(chrom, start, bin)
		if _, ok := binMap[curbin]; ok {
			binMap[curbin][1] = offset + length
		} else {
			bins = append(bins, curbin)
			binMap[curbin] = []int{offset, offset + length}
		}
		offset += length
	}
	for _, bin := range bins {
		fmt.Printf("%s\t%d\t%d\n", bin, binMap[bin][0], binMap[bin][1])
	}
	return nil
}

func RunTest(infile string) {
	IndexDatabse(infile, 1000)
}

func NewTestCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "test",
		Short: "Test",
		Run: func(cmd *cobra.Command, args []string) {
			input, _ := cmd.Flags().GetString("input")
			if input == "" {
				err := cmd.Help()
				if err != nil {
					log.Fatal(err)
				}
			} else {
				RunTest(input)
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "input")
	return cmd
}
