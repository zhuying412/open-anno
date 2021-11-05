package cmd

import (
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"grandanno/db"
	"grandanno/gene"
	"grandanno/gene_based/cnv"
	"grandanno/input"
	"grandanno/output"
	"log"
	"path"
)

func AnnoCnv(inputPath string, outputPath string) {
	// param
	db.InitDBParam()
	// chrom
	db.InitChromArray()
	// entrez_id
	db.InitSymbolToEntrezId()
	// refgene index
	refIndexes := db.GetRefIndexes()
	// refgene
	refgeneFile := path.Join(viper.GetString("db.path"), viper.GetString("db.refgene"))
	log.Printf("read %s\n", refgeneFile)
	refgeneMap := gene.ReadRefgeneFile(refgeneFile)
	// input
	log.Printf("read %s\n", inputPath)
	cnvs, infoMap := input.ReadCnvInputFile(inputPath)
	geneAnnoMap := cnv.NewGeneAnnoMap(cnvs, refgeneMap, refIndexes)
	log.Printf("write output %s\n", outputPath)
	outputs := make(output.CnvOutputs, 0)
	for _, c := range cnvs {
		outputs = append(outputs, output.CnvOutput{
			Cnv:       c,
			GeneAnnos: geneAnnoMap[c.SN()],
			OtherInfo: infoMap[c.SN()],
		})
	}
	output.CreateOutputFile(outputs, outputPath)
}

func NewCnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "cnv",
		Short: "Annotate CNV",
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
				AnnoCnv(cmd.Flag("input").Value.String(), cmd.Flag("output").Value.String())
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "SNV input file")
	cmd.Flags().StringP("output", "o", "", "Annotated ouput file")
	cmd.Flags().StringP("database", "d", "", "Database Path")
	return cmd
}
