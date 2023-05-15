package clinvar

import (
	"log"

	"fmt"
	"open-anno/pkg"
	"regexp"
	"strings"

	"github.com/brentp/bix"
	"github.com/brentp/faidx"
	"github.com/spf13/cobra"
)

type PrePathogenicMTParam struct {
	PrePathogenicParam
}

func (this PrePathogenicMTParam) GetVariants() ([]*pkg.SNV, error) {
	snvs, err := this.PrePathogenicParam.GetVariants()
	if err != nil {
		return []*pkg.SNV{}, err
	}
	MTSnvs := make([]*pkg.SNV, 0)
	for _, snv := range snvs {
		if snv.Chromosome == "chrM" {
			MTSnvs = append(MTSnvs, snv)
		}
	}
	return MTSnvs, nil
}

func (this PrePathogenicMTParam) Run() error {
	log.Printf("Read Gene: %s ...", this.Gene)
	err := pkg.InitGeneSymbolToID(this.Gene)
	if err != nil {
		return err
	}
	log.Printf("Open TABIX Handle ...")
	gpeTbx, err := bix.New(this.GenePred)
	if err != nil {
		return err
	}
	defer gpeTbx.Close()
	// 打开Genome
	genome, err := faidx.New(this.Genome)
	if err != nil {
		return err
	}
	snvs, err := this.GetVariants()
	if err != nil {
		return err
	}
	annoResult, err := this.RunAnno(snvs, gpeTbx, genome)
	if err != nil {
		return err
	}
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprint(writer, "Chrom\tPos\tRef\tAlt\tGene\tGeneID\tTranscript\tAAPos\tAARef\tAAAlt\n")
	for _, snv := range snvs {
		if info, ok := annoResult[snv.PK()]; ok {
			genes := strings.Split(info["GENE"].(string), ",")
			geneIds := strings.Split(info["GENE_ID"].(string), ",")
			geneDetails := strings.Split(info["DETAIL"].(string), ",")
			for i, gene := range genes {
				if gene != "." {
					for _, detail := range strings.Split(geneDetails[i], "|") {
						if detail != "." {
							info := strings.Split(detail, ":")
							if len(info) >= 5 {
								trans := info[1]
								re := regexp.MustCompile(`p\.([A-Z\*])(\d+)([A-Z\*])`)
								match := re.FindStringSubmatch(info[4])
								fmt.Println(info[4], match)
								fmt.Fprintf(writer, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", snv.Chromosome, snv.Pos, snv.Ref(), snv.Alt()[0], gene, geneIds[i], trans, match[2], match[1], match[3])
							}
						} else {
							fmt.Fprintf(writer, "%s\t%d\t%s\t%s\t%s\t%s\n", snv.Chromosome, snv.Pos, snv.Ref(), snv.Alt()[0], gene, geneIds[i])
						}
					}
				}
			}
		}
	}
	return err
}

func NewPrePathogenicMTCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "pathogenic_mt",
		Short: "Prepare MT Pathogenic Info from NCBI Clinvar VCF",
		Run: func(cmd *cobra.Command, args []string) {
			var param PrePathogenicMTParam
			param.Input, _ = cmd.Flags().GetString("input")
			param.Output, _ = cmd.Flags().GetString("output")
			param.GenePred, _ = cmd.Flags().GetString("genepred")
			param.Genome, _ = cmd.Flags().GetString("genome")
			param.Gene, _ = cmd.Flags().GetString("gene")
			param.AAshort, _ = cmd.Flags().GetBool("aashort")
			param.Exon, _ = cmd.Flags().GetBool("exon")
			err := param.Valid()
			if err != nil {
				cmd.Help()
				log.Fatal(err)
			}
			err = param.Run()
			if err != nil {
				log.Fatal(err)
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "Input Clinvar VCF File")
	cmd.Flags().StringP("output", "o", "", "Output File")
	cmd.Flags().StringP("genepred", "d", "", "Input GenePred File")
	cmd.Flags().StringP("genome", "G", "", "Input Genome Fasta File")
	cmd.Flags().StringP("gene", "g", "", "Input Gene Symbol To ID File")
	cmd.Flags().BoolP("aashort", "a", false, "Parameter Is AA Short")
	cmd.Flags().BoolP("exon", "e", false, "Parameter Is Exon")
	return cmd
}
