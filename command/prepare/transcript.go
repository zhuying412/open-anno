package prepare

import (
	"OpenAnno/command/viper"
	"OpenAnno/db/chromosome"
	"OpenAnno/db/transcript"
	"OpenAnno/db/transcript/index"
	"github.com/spf13/cobra"
	"log"
)

func NewPrepareTranscriptCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "transcript",
		Short: "Prepare required transcript files",
		Run: func(cmd *cobra.Command, args []string) {
			reference, _ := cmd.Flags().GetString("reference")
			refgene, _ := cmd.Flags().GetString("refgene")
			outdir, _ := cmd.Flags().GetString("outdir")
			database, _ := cmd.Flags().GetString("database")
			stream, _ := cmd.Flags().GetInt("stream")
			step, _ := cmd.Flags().GetInt("step")
			if reference == "" || refgene == "" || outdir == "" || database == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				viper.InitViper(database)
				chromosome.Init()
				transcript.Generate(refgene, reference, outdir, stream)
				index.Generate(refgene, outdir, stream, step)
			}
		},
	}
	cmd.Flags().StringP("reference", "r", "", "Reference Fasta File")
	cmd.Flags().StringP("refgene", "g", "", "RefGene File")
	cmd.Flags().StringP("outdir", "o", "", "Output Directory of Transcript Database")
	cmd.Flags().StringP("database", "d", "", "Database Path")
	cmd.Flags().IntP("stream", "l1", 3000, "Up/Down Stream Length")
	cmd.Flags().IntP("step", "l2", 300000, "Transcript Index Step Length")
	return cmd
}
