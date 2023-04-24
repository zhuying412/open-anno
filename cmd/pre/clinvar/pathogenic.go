package clinvar

// import (
// 	"log"

// 	"fmt"
// 	"io"
// 	"log"
// 	"open-anno/cmd/anno"
// 	"open-anno/pkg"
// 	"os"
// 	"path"
// 	"regexp"
// 	"strings"

// 	"github.com/brentp/vcfgo"
// 	"github.com/go-playground/validator/v10"
// 	"github.com/spf13/cobra"
// )

// type PrePathogenicParam struct {
// 	Input         string `validate:"required,pathexists"`
// 	Output        string `validate:"required"`
// 	GenePred      string `validate:"required,pathexists"`
// 	GenePredIndex string `validate:"required,pathexists"`
// 	Genome        string `validate:"required,pathexists"`
// 	Gene          string `validate:"required,pathexists"`
// 	GenomeIndex   string `validate:"required,pathexists"`
// 	AAshort       bool
// 	Exon          bool
// }

// func (this PrePathogenicParam) Valid() error {
// 	this.GenePredIndex = this.GenePred + ".tbi"
// 	this.GenomeIndex = this.Genome + ".fai"
// 	validate := validator.New()
// 	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
// 	err := validate.Struct(this)
// 	if err != nil {
// 		return err
// 	}
// 	outdir := path.Dir(this.Output)
// 	return os.MkdirAll(outdir, 0666)
// }

// func (this PrePathogenicParam) NewVcfWriter(writer io.WriteCloser, vcf string, infoKeys []string) (*vcfgo.Writer, error) {
// 	reader, err := pkg.NewIOReader(vcf)
// 	if err != nil {
// 		return &vcfgo.Writer{}, err
// 	}
// 	defer reader.Close()
// 	vcfReader, err := vcfgo.NewReader(reader, false)
// 	if err != nil {
// 		return &vcfgo.Writer{}, err
// 	}
// 	defer vcfReader.Close()
// 	vcfHeader := *vcfReader.Header
// 	vcfHeader.Infos = map[string]*vcfgo.Info{}
// 	return vcfgo.NewWriter(writer, &vcfHeader)
// }

// func (this PrePathogenicParam) FilterPathogenic(filteredFile string) error {
// 	reader, err := pkg.NewIOReader(this.Input)
// 	if err != nil {
// 		return err
// 	}
// 	defer reader.Close()
// 	vcfReader, err := vcfgo.NewReader(reader, false)
// 	if err != nil {
// 		return err
// 	}
// 	defer vcfReader.Close()
// 	writer, err := pkg.NewIOWriter(filteredFile)
// 	if err != nil {
// 		return err
// 	}
// 	defer writer.Close()
// 	vcfHeader := *vcfReader.Header
// 	vcfHeader.Infos = map[string]*vcfgo.Info{}
// 	vcfWriter, err := vcfgo.NewWriter(writer, vcfReader.Header)
// 	if err != nil {
// 		return err
// 	}
// 	for row := vcfReader.Read(); row != nil; row = vcfReader.Read() {
// 		row.Chromosome = "chr" + strings.ReplaceAll(row.Chromosome, "MT", "M")
// 		if len(row.Ref()) == 1 && len(row.Alt()[0]) == 1 {
// 			val, err := row.Info().Get("CLNREVSTAT")
// 			if err != nil {
// 				return err
// 			}
// 			star := 0
// 			revstat, ok := val.(string)
// 			if ok {
// 				revstat = strings.ToLower(revstat)
// 			} else {
// 				revstat = strings.ToLower(strings.Join(val.([]string), ";"))
// 			}
// 			if strings.Contains(revstat, "practice_guideline") {
// 				star = 4
// 			} else if strings.Contains(revstat, "reviewed_by_expert_panel") {
// 				star = 3
// 			} else if strings.Contains(revstat, "criteria_provided") {
// 				if strings.Contains(revstat, "multiple_submitters") {
// 					star = 2
// 				} else {
// 					star = 1
// 				}
// 			}
// 			if star >= 2 {
// 				val, err = row.Info().Get("CLNSIG")
// 				if err != nil {
// 					return err
// 				}
// 				clnsig, ok := val.(string)
// 				if ok {
// 					clnsig = strings.ToLower(clnsig)
// 				} else {
// 					clnsig = strings.ToLower(strings.Join(val.([]string), ";"))
// 				}
// 				if strings.Contains(clnsig, "pathogenic") && !strings.Contains(clnsig, "conflict") {
// 					row.Info_ = &vcfgo.InfoByte{}
// 					vcfWriter.WriteVariant(row)
// 				}
// 			}
// 		}
// 	}
// 	return nil
// }

// func (this PrePathogenicParam) AnnoGeneBase(filteredFile string) (string, error) {
// 	outFile := path.Join(path.Dir(this.Output), "clinvar.filtered.annotated.vcf")
// 	anno := anno.AnnoSnvParam{
// 		Input:       filteredFile,
// 		GenePred:    this.GenePred,
// 		Genome:      this.Genome,
// 		Gene:        this.Gene,
// 		Output:      outFile,
// 		AAshort:     true,
// 		Exon:        true,
// 		Concurrency: 16,
// 	}
// 	err := anno.Valid()
// 	if err != nil {
// 		return outFile, err
// 	}
// 	err = anno.Run()
// 	if err != nil {
// 		return outFile, err
// 	}
// 	return outFile, nil
// }

// func (this PrePathogenicParam) Run() error {
// 	filteredFile := path.Join(path.Dir(this.Output), "clinvar.filtered.vcf")
// 	if _, err := os.Stat(filteredFile); os.IsNotExist(err) {
// 		err := this.FilterPathogenic(filteredFile)
// 		if err != nil {
// 			return err
// 		}
// 	}
// 	annoResult, err := this.AnnoGeneBase(filteredFile)
// 	if err != nil {
// 		return err
// 	}
// 	reader, err := pkg.NewIOReader(this.Input)
// 	if err != nil {
// 		return err
// 	}
// 	defer reader.Close()
// 	vcfReader, err := vcfgo.NewReader(reader, false)
// 	if err != nil {
// 		return err
// 	}
// 	defer vcfReader.Close()
// 	writer, err := pkg.NewIOWriter(this.Output)
// 	if err != nil {
// 		return err
// 	}
// 	defer writer.Close()
// 	fmt.Fprint(writer, "Chrom\tPos\tRef\tAlt\tGene\tGeneID\tTranscript\tAAPos\tAARef\tAAAlt\n")
// 	for row := vcfReader.Read(); row != nil; row = vcfReader.Read() {
// 		row.Chromosome = "chr" + row.Chromosome
// 		snv := &pkg.Variant{Variant: *row}
// 		igeneinfo, err := row.Info().Get("GENEINFO")
// 		if err != nil {
// 			continue
// 		}
// 		geneinfos := strings.Split(strings.Join(pkg.Interface2Array[string](igeneinfo), "|"), "|")
// 		clnGeneIds := make([]string, len(geneinfos))
// 		for i, geneinfo := range geneinfos {
// 			clnGeneIds[i] = strings.Split(geneinfo, ":")[1]
// 		}
// 		if annoInfos, ok := annoResult.AnnoInfos[snv.PK()]; ok && len(annoInfos) > 0 {
// 			genes := strings.Split(annoInfos["GENE"].(string), ",")
// 			geneIds := strings.Split(annoInfos["GENE_ID"].(string), ",")
// 			geneDetails := strings.Split(annoInfos["DETAIL"].(string), ",")
// 			for i, gene := range genes {
// 				geneId := geneIds[i]
// 				if pkg.FindArr(clnGeneIds, geneId) != -1 {
// 					for _, detail := range strings.Split(geneDetails[i], "|") {
// 						info := strings.Split(detail, ":")
// 						if len(info) >= 5 {
// 							trans := info[1]
// 							re := regexp.MustCompile(`p\.([A-Z\*])(\d+)([A-Z\*])`)
// 							match := re.FindStringSubmatch(info[4])
// 							fmt.Println(info[4], match)
// 							fmt.Fprintf(writer, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", row.Chromosome, row.Pos, row.Ref(), row.Alt()[0], gene, geneIds[i], trans, match[2], match[1], match[3])
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
// 	return err
// }

// func NewPrePathogenicCmd() *cobra.Command {
// 	cmd := &cobra.Command{
// 		Use:   "pathogenic",
// 		Short: "Prepare Pathogenic Info from NCBI Clinvar VCF",
// 		Run: func(cmd *cobra.Command, args []string) {
// 			var param PrePathogenicParam
// 			param.Input, _ = cmd.Flags().GetString("input")
// 			param.Output, _ = cmd.Flags().GetString("output")
// 			param.GenePred, _ = cmd.Flags().GetString("genepred")
// 			param.Genome, _ = cmd.Flags().GetString("genome")
// 			param.Gene, _ = cmd.Flags().GetString("gene")
// 			param.AAshort, _ = cmd.Flags().GetBool("aashort")
// 			param.Exon, _ = cmd.Flags().GetBool("exon")
// 			err := param.Valid()
// 			if err != nil {
// 				cmd.Help()
// 				log.Fatal(err)
// 			}
// 			err = param.Run()
// 			if err != nil {
// 				log.Fatal(err)
// 			}
// 		},
// 	}
// 	cmd.Flags().StringP("input", "i", "", "Input Clinvar VCF File")
// 	cmd.Flags().StringP("output", "o", "", "Output File")
// 	cmd.Flags().StringP("genepred", "d", "", "Input GenePred File")
// 	cmd.Flags().StringP("genome", "G", "", "Input Genome Fasta File")
// 	cmd.Flags().StringP("gene", "g", "", "Input Gene Symbol To ID File")
// 	cmd.Flags().BoolP("aashort", "a", false, "Parameter Is AA Short")
// 	cmd.Flags().BoolP("exon", "e", false, "Parameter Is Exon")
// 	return cmd
// }
