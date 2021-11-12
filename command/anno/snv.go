package anno

import (
	"OpenAnno/anno"
	"OpenAnno/anno/gene/snv"
	viper2 "OpenAnno/command/viper"
	"OpenAnno/db/chromosome"
	"OpenAnno/db/gene"
	"OpenAnno/variant"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"log"
)

func runSnvGeneAnno(annoMap *map[string]map[string]anno.IAnno, snvs variant.Snvs, databaseDir string, chrom string) {
	for key, val := range viper.GetStringMapString("transcript") {
		transMap, transIndexes := readTranscriptDir(databaseDir, val, chrom)
		log.Printf("run annotate of %s", chrom)
		geneAnnoMap := snv.RunAnnotate(snvs, transMap, transIndexes)
		for _, _snv := range snvs {
			(*annoMap)[_snv.SN()][key] = geneAnnoMap[_snv.SN()]
		}
	}
}

func AnnoSnv(inputFile string, outputFile string) {
	databaseDir := viper.GetString("db.path")
	// snv
	log.Printf("read %s", inputFile)
	snvs, infoMap := variant.ReadSnvFile(inputFile)
	annoMap := make(map[string]map[string]anno.IAnno)
	for _, chrom := range chromosome.ChromList {
		subSnvs := snvs.FilterByChrom(chrom.Name)
		runSnvGeneAnno(&annoMap, subSnvs, databaseDir, chrom.Name)
		runFilterAnno(&annoMap, subSnvs, databaseDir, chrom.Name)
		runRegionAnno(&annoMap, subSnvs, databaseDir, chrom.Name)
	}
	writeResult(snvs, annoMap, infoMap, outputFile)
}

func NewAnnoSnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "snv",
		Short: "Annotate SNV",
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
				gene.Init()
				snv.Init()
				AnnoSnv(input, output)
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "SNV input file")
	cmd.Flags().StringP("output", "o", "", "Annotated ouput file")
	cmd.Flags().StringP("database", "d", "", "Database Path")
	return cmd
}
