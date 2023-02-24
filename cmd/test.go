package cmd

import (
	"fmt"
	"log"
	"open-anno/pkg"

	"github.com/brentp/vcfgo"
	"github.com/spf13/cobra"
)

func RunTest(infile string) {
	reader, err := pkg.NewIOReader(infile)
	if err != nil {
		log.Fatalln(err)
	}
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		log.Fatalln(err)
	}
	// fmt.Println(vcfReader.AddInfoToHeader())
	for variant := vcfReader.Read(); variant != nil; {
		fmt.Println(variant.)
	}
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
