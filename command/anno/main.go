package anno

import (
	"OpenAnno/anno"
	"OpenAnno/anno/database"
	"OpenAnno/db/transcript"
	"OpenAnno/db/transcript/index"
	"OpenAnno/variant"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"log"
	"path"
)

func readTranscriptDir(databaseDir string, name string, chrom string) (transcript.TranscriptMap, index.TranscriptIndexes) {
	transFile := path.Join(databaseDir, name, "chr"+chrom+".json")
	log.Printf("read %s", transFile)
	transMap := transcript.ReadTranscriptJSON(transFile)
	transIndexFile := path.Join(databaseDir, name, "chr"+chrom+".idx.json")
	log.Printf("read %s", transIndexFile)
	transIndexes := index.ReadTranscriptIndexJSON(transIndexFile)
	return transMap, transIndexes
}

func runFilterAnno(annoMap *map[string]map[string]anno.IAnno, variants variant.IVariants, databaseDir string, chrom string) {
	for key, val := range viper.GetStringMapString("filter_based") {
		databaseFile := path.Join(databaseDir, val, "chr"+chrom+".txt")
		filterAnnoMap := database.RunFilterAnnotate(variants, databaseFile)
		for i := 0; i < variants.Len(); i++ {
			sn := variants.GetVariant(i).SN()
			(*annoMap)[sn][key] = filterAnnoMap[sn]
		}
	}
}

func runRegionAnno(annoMap *map[string]map[string]anno.IAnno, variants variant.IVariants, databaseDir string, chrom string) {
	for key, val := range viper.GetStringMapString("region_based") {
		databaseFile := path.Join(databaseDir, val, "chr"+chrom+".txt")
		regionAnnoMap := database.RunRegionAnnotate(variants, databaseFile)
		for i := 0; i < variants.Len(); i++ {
			sn := variants.GetVariant(i).SN()
			(*annoMap)[sn][key] = regionAnnoMap[sn]
		}
	}

}

func writeResult(variants variant.IVariants, annoMap map[string]map[string]anno.IAnno, infoMap variant.OtherInfoMap, outputFile string) {
	log.Printf("write result to %s", outputFile)
	results := make([]anno.Result, 0)
	for i := 0; i < variants.Len(); i++ {
		_variant := variants.GetVariant(i)
		results = append(results, anno.Result{
			Variant:    _variant,
			Annotation: annoMap[_variant.SN()],
			OtherInfo:  infoMap[_variant.SN()],
		})
	}
	anno.CreateResultJSON(results, outputFile)
}

func NewAnnoCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "anno",
		Short: "Run annotate",
	}
	cmd.AddCommand(NewAnnoSnvCmd())
	cmd.AddCommand(NewAnnoCnvCmd())
	return cmd
}
