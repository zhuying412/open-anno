package clinvar

import (
	"log"
	"sync"

	"fmt"
	"open-anno/anno"
	"open-anno/pkg"
	"os"
	"path"
	"regexp"
	"strings"

	"github.com/brentp/bix"
	"github.com/brentp/faidx"
	"github.com/brentp/vcfgo"
	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
	"github.com/syndtr/goleveldb/leveldb"
)

type PrePathogenicParam struct {
	Input         string `validate:"required,pathexists"`
	Output        string `validate:"required"`
	GenePred      string `validate:"required,pathexists"`
	GenePredIndex string `validate:"required,pathexists"`
	Genome        string `validate:"required,pathexists"`
	Gene          string `validate:"required,pathexists"`
	GenomeIndex   string `validate:"required,pathexists"`
	AAshort       bool
	Exon          bool
}

func (this PrePathogenicParam) Valid() error {
	this.GenePredIndex = this.GenePred + ".tbi"
	this.GenomeIndex = this.Genome + ".fai"
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this PrePathogenicParam) GetVariants() ([]*pkg.SNV, error) {
	snvs := make([]*pkg.SNV, 0)
	reader, err := pkg.NewIOReader(this.Input)
	if err != nil {
		return []*pkg.SNV{}, err
	}
	defer reader.Close()
	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return []*pkg.SNV{}, err
	}
	defer vcfReader.Close()
	for row := vcfReader.Read(); row != nil; row = vcfReader.Read() {
		snv := &pkg.SNV{*row}
		snv.Chromosome = "chr" + strings.ReplaceAll(snv.Chromosome, "MT", "M")
		if len(snv.Ref()) == 1 && len(snv.Alt()[0]) == 1 {
			val, err := snv.Info().Get("CLNREVSTAT")
			if err != nil {
				return []*pkg.SNV{}, err
			}
			star := 0
			revstat, ok := val.(string)
			if ok {
				revstat = strings.ToLower(revstat)
			} else {
				revstat = strings.ToLower(strings.Join(val.([]string), ";"))
			}
			if strings.Contains(revstat, "practice_guideline") {
				star = 4
			} else if strings.Contains(revstat, "reviewed_by_expert_panel") {
				star = 3
			} else if strings.Contains(revstat, "criteria_provided") {
				if strings.Contains(revstat, "multiple_submitters") {
					star = 2
				} else {
					star = 1
				}
			}
			if star >= 2 {
				val, err = snv.Info().Get("CLNSIG")
				if err != nil {
					return []*pkg.SNV{}, err
				}
				clnsig, ok := val.(string)
				if ok {
					clnsig = strings.ToLower(clnsig)
				} else {
					clnsig = strings.ToLower(strings.Join(val.([]string), ";"))
				}
				if strings.Contains(clnsig, "pathogenic") && !strings.Contains(clnsig, "conflict") {
					snvs = append(snvs, snv)
				}
			}
		}
	}
	return snvs, nil
}

func (this PrePathogenicParam) RunAnno(snvs []*pkg.SNV, gpeTbx *bix.Bix, genome *faidx.Faidx) (map[string]map[string]any, error) {
	snvChan := make(chan *pkg.SNV, len(snvs))
	for _, snv := range snvs {
		snvChan <- snv
	}
	close(snvChan)
	var wg sync.WaitGroup
	resChan := make(chan anno.AnnoInfo, len(snvs))
	for i := 0; i <= 80; i++ {
		wg.Add(1)
		go anno.AnnoSnvWorker(snvChan, gpeTbx, []*bix.Bix{}, []*bix.Bix{}, []string{}, []*leveldb.DB{}, genome, 0.7, resChan, &wg)
	}
	go func() {
		wg.Wait()
		close(resChan)
	}()
	results := make(map[string]map[string]any)
	for res := range resChan {
		if res.Error != nil {
			return results, res.Error
		}
		results[res.PK] = res.Data
	}
	return results, nil
}

func (this PrePathogenicParam) Run() error {
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
		igeneinfo, err := snv.Info().Get("GENEINFO")
		if err != nil {
			continue
		}
		geneinfos := strings.Split(strings.Join(pkg.Interface2Array[string](igeneinfo), "|"), "|")
		clnGeneIds := make([]string, len(geneinfos))
		for i, geneinfo := range geneinfos {
			clnGeneIds[i] = strings.Split(geneinfo, ":")[1]
		}
		if info, ok := annoResult[snv.PK()]; ok {
			genes := strings.Split(info["GENE"].(string), ",")
			geneIds := strings.Split(info["GENE_ID"].(string), ",")
			geneDetails := strings.Split(info["DETAIL"].(string), ",")
			for i, gene := range genes {
				if gene != "." {
					geneId := geneIds[i]
					if pkg.FindArr(clnGeneIds, geneId) != -1 {
						for _, detail := range strings.Split(geneDetails[i], "|") {
							info := strings.Split(detail, ":")
							if len(info) >= 5 {
								trans := info[1]
								re := regexp.MustCompile(`p\.([A-Z\*])(\d+)([A-Z\*])`)
								match := re.FindStringSubmatch(info[4])
								fmt.Println(info[4], match)
								fmt.Fprintf(writer, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", snv.Chromosome, snv.Pos, snv.Ref(), snv.Alt()[0], gene, geneIds[i], trans, match[2], match[1], match[3])
							}
						}
					}
				}
			}
		}
	}
	return err
}

func NewPrePathogenicCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "pathogenic",
		Short: "Prepare Pathogenic Info from NCBI Clinvar VCF",
		Run: func(cmd *cobra.Command, args []string) {
			var param PrePathogenicParam
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
