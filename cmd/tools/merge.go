package tools

import (
	"log"
	"open-anno/tools"

	"github.com/spf13/cobra"
)

func NewMergeCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "merge",
		Short: "Merge annotations",
		Run: func(cmd *cobra.Command, args []string) {
			genebased, _ := cmd.Flags().GetString("genebased")
			otherbaseds, _ := cmd.Flags().GetStringArray("otherbased")
			outfile, _ := cmd.Flags().GetString("outfile")
			if genebased == "" || outfile == "" || len(otherbaseds) == 0 {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				tools.MergeAnno(outfile, genebased, otherbaseds...)
			}
		},
	}
	cmd.Flags().StringP("genebased", "g", "", "Genebased Annotation File")
	cmd.Flags().StringArrayP("otherbased", "d", []string{}, "Filterbased or Regionbased Annotation File")
	cmd.Flags().StringP("outfile", "o", "", "Output File")
	return cmd
}
