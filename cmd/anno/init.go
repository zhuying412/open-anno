package anno

import (
	"fmt"
	"log"
	"open-anno/pkg/io"
	"os"
	"strings"

	"github.com/spf13/cobra"
)

func NewAnnoCmd(varType string) *cobra.Command {
	varType = strings.ToLower(varType)
	cmd := &cobra.Command{
		Use:   varType,
		Short: fmt.Sprintf("Annotate for %s", strings.ToUpper(varType)),
		Run: func(cmd *cobra.Command, args []string) {
			avinput, _ := cmd.Flags().GetString("avinput")
			dbpath, _ := cmd.Flags().GetString("dbpath")
			dbnames, _ := cmd.Flags().GetString("dbnames")
			dbtypes, _ := cmd.Flags().GetString("dbtypes")
			builder, _ := cmd.Flags().GetString("builder")
			outprefix, _ := cmd.Flags().GetString("outprefix")
			if avinput == "" || dbpath == "" || dbnames == "" || dbtypes == "" || builder == "" || outprefix == "" {
				fmt.Println(avinput, dbpath, dbnames, dbtypes, builder, outprefix)
				err := cmd.Help()
				if err != nil {
					log.Panic(err)
				}
			} else {
				if varType != "snv" && varType != "cnv" {
					log.Fatalln("only 'snv' or 'cnv'")
				}
				dbNames := strings.Split(dbnames, ",")
				dbTypes := strings.Split(dbtypes, ",")
				var genebasedOut string
				var otherbasedOuts []string
				errChan := make(chan error, len(dbNames))
				for i := 0; i < len(dbNames); i++ {
					dbName := strings.TrimSpace(dbNames[i])
					dbType := strings.TrimSpace(dbTypes[i])
					outfile := fmt.Sprintf("%s.%s.anno.txt", outprefix, dbName)
					if dbType == "g" {
						if varType == "snv" {
							aashort, _ := cmd.Flags().GetBool("aashort")
							go RunAnnoSnvGeneBased(avinput, dbpath, dbName, builder, outfile, aashort, errChan)
						}
						if varType == "cnv" {
							go RunAnnoCnvGeneBased(avinput, dbpath, dbName, builder, outfile, errChan)
						}
						genebasedOut = outfile
					}
					if dbType == "f" {
						go RunAnnoFilterBased(avinput, dbpath, dbName, builder, outfile, errChan)
						otherbasedOuts = append(otherbasedOuts, outfile)
					}
					if dbType == "r" {
						overlap, _ := cmd.Flags().GetFloat64("overlap")
						go RunAnnoRegionBased(avinput, dbpath, dbName, builder, overlap, outfile, errChan)
						otherbasedOuts = append(otherbasedOuts, outfile)
					}
				}
				for i := 0; i < len(dbNames); i++ {
					err := <-errChan
					if err != nil {
						log.Fatal(err)
					}
				}
				outfile := fmt.Sprintf("%s.anno.txt", outprefix)
				io.MergeAnno(outfile, avinput, genebasedOut, otherbasedOuts...)
				if isClean, _ := cmd.Flags().GetBool("clean"); isClean {
					os.Remove(genebasedOut)
					for _, otherbasedOut := range otherbasedOuts {
						os.Remove(otherbasedOut)
					}
				}
			}
		},
	}
	cmd.Flags().StringP("avinput", "i", "", "Annotated Variants Input File")
	cmd.Flags().StringP("outprefix", "o", "", "Output Prefix")
	cmd.Flags().StringP("dbpath", "d", "", "Database Directory")
	cmd.Flags().StringP("dbnames", "n", "", "Database Names")
	cmd.Flags().StringP("dbtypes", "t", "", "Database Types")
	cmd.Flags().StringP("builder", "b", "hg19", "Database Builder")
	cmd.Flags().BoolP("clean", "c", false, "Database Builder")
	if varType == "snv" {
		cmd.Flags().BoolP("aashort", "s", false, "Database Builder")
	}
	if varType == "cnv" {
		cmd.Flags().Float64P("overlap", "p", 0.75, "CNV Overlap Threshold")
	}
	return cmd
}
