package cmd

import (
	"log"

	"github.com/spf13/cobra"
)

func RunTest(infile string) {
	return
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
