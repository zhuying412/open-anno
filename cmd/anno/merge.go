package anno

import (
	"log"
	"open-anno/pkg/io"

	"github.com/spf13/cobra"
)

func NewMergeCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "merge",
		Short: "Merge Annotation",
		Run: func(cmd *cobra.Command, args []string) {
			avinput, _ := cmd.Flags().GetString("avinput")
			genebased, _ := cmd.Flags().GetString("genebased")
			otherbaseds, _ := cmd.Flags().GetStringArray("otherbaseds")
			outfile, _ := cmd.Flags().GetString("outfile")
			if avinput == "" || genebased == "" || outfile == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				io.MergeAnno(outfile, avinput, genebased, otherbaseds...)
			}
		},
	}
	cmd.Flags().StringP("avinput", "i", "", "Annotated Variants Input File")
	cmd.Flags().StringP("genebased", "g", "", "GeneBased Annotation File")
	cmd.Flags().StringArrayP("otherbaseds", "d", []string{}, "FilterBased or RegionBased Annotation Files")
	cmd.Flags().StringP("outfile", "o", "", "Output Merged Annotation File")
	return cmd
}
