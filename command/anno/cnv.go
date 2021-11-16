package anno

import (
	"OpenAnno/anno"
	"OpenAnno/anno/gene/cnv"
	viper2 "OpenAnno/command/viper"
	"OpenAnno/db/chromosome"
	"OpenAnno/variant"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"log"
)

func AnnoCnv(inputFile string, outputFile string) {
	databaseDir := viper.GetString("db.path")
	// snv
	log.Printf("read %s", inputFile)
	cnvs, infoMap := variant.ReadCnvFile(inputFile)
	annoMap := make(map[string]map[string]anno.IAnno)
	for _, chrom := range chromosome.ChromList {
		subCnvs := cnvs.FilterByChrom(chrom.Name)
		for key, val := range viper.GetStringMapString("transcript") {
			transMap, transIndexes := readTranscriptDir(databaseDir, val, chrom.Name)
			log.Printf("run annotate of %s", chrom.Name)
			geneAnnoMap := cnv.RunAnnotate(subCnvs, transMap, transIndexes)
			for _, _cnv := range subCnvs {
				annoMap[_cnv.SN()][key] = geneAnnoMap[_cnv.SN()]
			}
		}
		runRegionAnno(&annoMap, subCnvs, databaseDir, chrom.Name)
	}
	writeResult(cnvs, annoMap, infoMap, outputFile)
}

func NewAnnoCnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "cnv",
		Short: "Annotate CNV",
		Run: func(cmd *cobra.Command, args []string) {
			input, _ := cmd.Flags().GetString("input")
			output, _ := cmd.Flags().GetString("output")
			database, _ := cmd.Flags().GetString("database")

			if input == "" || output == "" || database == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				viper2.InitViper(database)
				chromosome.Init()
				AnnoCnv(input, output)
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "SNV input file")
	cmd.Flags().StringP("output", "o", "", "Annotated ouput file")
	cmd.Flags().StringP("database", "d", "", "Database Path")
	return cmd
}
