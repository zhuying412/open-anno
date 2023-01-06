package tools

import (
	"fmt"
	"log"
	"open-anno/pkg"
	"os"
	"path"
	"sort"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
)

type GeneAnno struct {
	Snv     string
	Gene    string
	GeneID  string
	Regions []string
	Events  []string
	Details []string
}

func (this *GeneAnno) AddEvent(event string) {
	if event != "." && pkg.FindArr(this.Events, event) < 0 {
		this.Events = append(this.Events, event)
	}
}
func (this *GeneAnno) AddRegion(region string) {
	if region != "." && pkg.FindArr(this.Regions, region) < 0 {
		this.Regions = append(this.Regions, region)
	}
}
func (this *GeneAnno) AddDetail(detail string) {
	if detail != "." && pkg.FindArr(this.Details, detail) < 0 {
		this.Details = append(this.Details, detail)
	}
}

func (this GeneAnno) Region() string {
	var regions1, regions2 []string
	for _, region := range this.Regions {
		switch region {
		case "exonic", "splicing", "exonic_splicing", "transcript":
			regions1 = append(regions1, region)
		case "ncRNA", "UTR3", "UTR5", "intronic":
			regions2 = append(regions2, region)
		}
	}
	if len(regions1) > 0 {
		return strings.Join(regions1, ",")
	}
	if len(regions2) > 0 {
		return strings.Join(regions2, ",")
	}
	return "."
}

func (this GeneAnno) Event() string {
	if len(this.Events) > 0 {
		return strings.Join(this.Events, ",")
	}
	return "."
}

func (this GeneAnno) Detail() string {
	if len(this.Details) > 0 {
		return strings.Join(this.Details, ",")
	}
	return "."
}

func (this GeneAnno) Text() string {
	return fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s", this.Snv, this.Gene, this.GeneID, this.Event(), this.Region(), this.Detail())
}

type AggsParam struct {
	Input  string `validate:"required,pathexists"`
	Output string `validate:"required"`
}

func (this AggsParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	outdir := path.Dir(this.Output)
	return os.MkdirAll(outdir, 0666)
}

func (this AggsParam) ReadTransAnno() (map[string]GeneAnno, string, error) {
	var header string
	geneAnnos := make(map[string]GeneAnno)
	reader, err := pkg.NewIOReader(this.Input)
	if err != nil {
		return geneAnnos, header, err
	}
	defer reader.Close()
	scanner := pkg.NewIOScanner(reader)
	scanner.Scan()
	header = scanner.Text()
	for scanner.Scan() {
		row := strings.Split(scanner.Text(), "\t")
		snv := strings.Join(row[0:5], "\t")
		gene, entrezId, event, region, detail := row[5], row[6], row[7], row[8], row[9]
		pk := fmt.Sprintf("%s\t%s", snv, gene)
		geneAnno, ok := geneAnnos[pk]
		if !ok {
			geneAnno = GeneAnno{Snv: snv, Gene: gene, GeneID: entrezId}
		}
		geneAnno.AddEvent(event)
		geneAnno.AddRegion(region)
		geneAnno.AddDetail(detail)
		geneAnnos[pk] = geneAnno
	}
	return geneAnnos, header, nil
}

func (this AggsParam) Run() error {
	geneAnnos, header, err := this.ReadTransAnno()
	if err != nil {
		return nil
	}
	geneAnnoTexts := make([]string, len(geneAnnos))
	i := 0
	for pk, geneAnno := range geneAnnos {
		fmt.Println(pk)
		geneAnnoTexts[i] = geneAnno.Text()
		i++
	}
	sort.Strings(geneAnnoTexts)
	writer, err := pkg.NewIOWriter(this.Output)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprintf(writer, "%s\n", header)
	for _, text := range geneAnnoTexts {
		fmt.Fprintf(writer, "%s\n", text)
	}
	return nil
}

func NewAggsCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "aggs",
		Short: "Aggregate Snv GeneBased Annotation",
		Run: func(cmd *cobra.Command, args []string) {
			var param AggsParam
			param.Input, _ = cmd.Flags().GetString("gbanno")
			param.Output, _ = cmd.Flags().GetString("output")
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
	cmd.Flags().StringP("gbanno", "i", "", "Input Annotate Result File of GeneBased")
	cmd.Flags().StringP("output", "o", "", "Output Aggregated Annotation File")
	return cmd
}
