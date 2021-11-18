package anno

import (
	"OpenAnno/anno"
	"OpenAnno/anno/database"
	transcript2 "OpenAnno/db/transcript"
	"OpenAnno/pkg/transcript"
	"OpenAnno/pkg/transcript/index"
	"OpenAnno/pkg/variant"
	"OpenAnno/run"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"log"
	"path"
)

func readTranscriptDir(databaseDir string, name string, chrom string) (transcript.TranscriptMap, index.TranscriptIndexes) {
	transFile := path.Join(databaseDir, name, "chr"+chrom+".json")
	transIndexFile := path.Join(databaseDir, name, "chr"+chrom+".idx.json")
	log.Printf("read gene anno file %s", transFile)
	transMap := transcript2.ReadTranscriptJSON(transFile)
	transIndexes := transcript2.ReadTranscriptIndexJSON(transIndexFile)
	return transMap, transIndexes
}

func runFilterAnno(annoMap *map[string]map[string]anno.IAnno, variants variant.IVariants, databaseDir string, chrom string) {
	for key, val := range viper.GetStringMapString("db.filter_based") {
		databaseFile := path.Join(databaseDir, val, "chr"+chrom+".txt")
		log.Printf("run filter based anno with %s", databaseFile)
		filterAnnoMap := database.RunFilterAnnotate(variants, databaseFile)
		for i := 0; i < variants.Len(); i++ {
			sn := variants.GetVariant(i).SN()
			(*annoMap)[sn][key] = filterAnnoMap[sn]
		}
	}
}

func runRegionAnno(annoMap *map[string]map[string]anno.IAnno, variants variant.IVariants, databaseDir string, chrom string) {
	for key, val := range viper.GetStringMapString("db.region_based") {
		databaseFile := path.Join(databaseDir, val, "chr"+chrom+".txt")
		log.Printf("run region based anno with %s", databaseFile)
		regionAnnoMap := database.RunRegionAnnotate(variants, databaseFile)
		for i := 0; i < variants.Len(); i++ {
			sn := variants.GetVariant(i).SN()
			(*annoMap)[sn][key] = regionAnnoMap[sn]
		}
	}

}

func writeResult(variants variant.IVariants, annoMap map[string]map[string]anno.IAnno, infoMap run.OtherInfoMap, outputFile string) {
	log.Printf("write result to %s", outputFile)
	results := make([]run.Result, 0)
	for i := 0; i < variants.Len(); i++ {
		_variant := variants.GetVariant(i)
		results = append(results, run.Result{
			Variant:    _variant,
			Annotation: annoMap[_variant.SN()],
			OtherInfo:  infoMap[_variant.SN()],
		})
	}
	run.CreateResultJSON(results, outputFile)
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
