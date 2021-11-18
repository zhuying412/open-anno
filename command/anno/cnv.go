package anno

import (
	"OpenAnno/anno"
	"OpenAnno/anno/gene/cnv"
	viper2 "OpenAnno/command/viper"
	"OpenAnno/db"
	"OpenAnno/run"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"log"
)

func AnnoCnv(inputFile string, outputFile string) {
	databaseDir := viper.GetString("db.path")
	// snv
	log.Printf("read %s", inputFile)
	cnvs, infoMap := run.ReadCnvFile(inputFile)
	annoMap := make(map[string]map[string]anno.IAnno)
	for _, chrom := range db.ChromList {
		subCnvs := cnvs.FilterByChrom(chrom.Name)
		for key, val := range viper.GetStringMapString("transcript") {
			transMap, transIndexes := readTranscriptDir(databaseDir, val, chrom.Name)
			log.Printf("run gene anno in chr%s", chrom)
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
			inputFile, _ := cmd.Flags().GetString("input")
			outputFile, _ := cmd.Flags().GetString("output")
			database, _ := cmd.Flags().GetString("database")

			if inputFile == "" || outputFile == "" || database == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				viper2.InitViper(database)
				db.InitChrom()
				AnnoCnv(inputFile, outputFile)
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "SNV input file")
	cmd.Flags().StringP("output", "o", "", "Annotated ouput file")
	cmd.Flags().StringP("database", "d", "", "Database Path")
	return cmd
}
