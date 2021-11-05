package cmd

import (
	"github.com/spf13/cobra"
	"grandanno/db"
	"grandanno/gene"
	"grandanno/gene_based/snv"
	"grandanno/input"
	"grandanno/output"
	"log"
)

func AnnoSnv(inputPath string, outputPath string) {
	// param
	db.InitDBParam()
	snv.InitSnvParam()
	// chrom
	db.InitChromArray()
	// entrez_id
	db.InitSymbolToEntrezId()
	// refgene index
	refIndexes := db.GetRefIndexes()
	// refgene
	refgeneMap := gene.GetRefgeneMap()
	// avinput
	log.Printf("read %s\n", inputPath)
	snvs, infoMap := input.ReadSnvInputFile(inputPath)
	geneAnnoMap := snv.NewGeneAnnoMap(snvs, refgeneMap, refIndexes)
	log.Printf("write output %s\n", outputPath)
	outputs := make(output.SnvOutputs, 0)
	for _, s := range snvs {
		outputs = append(outputs, output.SnvOutput{
			Snv:       s,
			GeneAnnos: geneAnnoMap[s.SN()],
			OtherInfo: infoMap[s.SN()],
		})
	}
	output.CreateOutputFile(outputs, outputPath)
}

func NewSnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "snv",
		Short: "Annotate SNV",
		Run: func(cmd *cobra.Command, args []string) {
			log.Println(cmd.Flag("database").Value.String())
			if cmd.Flag("database").Value.String() == "" ||
				cmd.Flag("input").Value.String() == "" ||
				cmd.Flag("output").Value.String() == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				InitViper(cmd.Flag("database").Value.String())
				AnnoSnv(cmd.Flag("input").Value.String(), cmd.Flag("output").Value.String())
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "SNV input file")
	cmd.Flags().StringP("output", "o", "", "Annotated ouput file")
	cmd.Flags().StringP("database", "d", "", "Database Path")
	return cmd
}
