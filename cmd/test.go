package cmd

import (
	"fmt"
	"log"

	"github.com/brentp/bix"
	"github.com/spf13/cobra"
)

type TestStruct struct {
	chrom string
	pos   uint32
}

func (this TestStruct) Chrom() string {
	return this.chrom
}

func (this TestStruct) Start() uint32 {
	return this.pos - 1
}
func (this TestStruct) End() uint32 {
	return this.pos
}

func RunTest(infile string) error {
	tbx, err := bix.New(infile)
	if err != nil {
		return err
	}
	defer tbx.Close()
	query, err := tbx.Query(TestStruct{"chr1", 11874})
	if err != nil {
		return err
	}
	for v, e := query.Next(); e == nil; v, e = query.Next() {
		fmt.Println(v.Chrom(), v.Start(), v.End())
	}
	return err
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
