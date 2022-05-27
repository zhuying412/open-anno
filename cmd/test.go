package cmd

import (
	"bufio"
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

func RunTest(infile string) {
	reader, _ := os.Open(infile)
	scanner := bufio.NewScanner(reader)
	var i int
	var curOffset int64
	for scanner.Scan() {
		i++
		if i%10000 == 0 {
			fmt.Println(i)
			if i%50000 == 0 {
				curOffset, _ = reader.Seek(0, os.SEEK_CUR)
				fmt.Println(curOffset)
				fmt.Println(scanner.Text())
				break
			}
		}
	}
	reader.Close()
	fmt.Println("=====================")
	reader, _ = os.Open(infile)
	// scanner = bufio.NewScanner(reader)
	reader.Seek(curOffset, os.SEEK_SET)
	// scanner.Scan()
	// fmt.Println(scanner.Text())
	bufReader := bufio.NewReader(reader)
	line, _, _ := bufReader.ReadLine()
	fmt.Printf("%s", line)
	reader.Close()
}

func NewTestCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "test",
		Short: "Test",
		Run: func(cmd *cobra.Command, args []string) {
			input, _ := cmd.Flags().GetStringArray("input")
			for _, i := range input {
				fmt.Println(i)
			}
			// if input == "" {
			// 	err := cmd.Help()
			// 	if err != nil {
			// 		log.Panic(err)
			// 	}
			// } else {
			// 	RunTest(input)
			// }
		},
	}
	cmd.Flags().StringArrayP("input", "i", []string{}, "input")
	return cmd
}
