package tools

import (
	"fmt"
	"io"
	"log"
	"open-anno/pkg"
	"os"
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

func (this *GeneAnno) AddInfo(event, region, detail string) {
	if event != "." && pkg.FindArr(this.Events, event) < 0 {
		this.Events = append(this.Events, event)
	}
	if region != "." && pkg.FindArr(this.Regions, region) < 0 {
		this.Regions = append(this.Regions, region)
	}
	if detail != "." && pkg.FindArr(this.Details, detail) < 0 {
		this.Details = append(this.Details, detail)
	}
}

func (this GeneAnno) Region() string {
	var regions1, regions2, regions3, regions4 []string
	for _, region := range this.Regions {
		switch region {
		case "exonic", "splicing", "exonic_splicing", "transcript":
			regions1 = append(regions1, region)
		case "ncRNA":
			regions2 = append(regions2, region)
		case "UTR3", "UTR5":
			regions3 = append(regions3, region)
		case "intronic":
			regions4 = append(regions4, region)
		}
	}
	if len(regions1) > 0 {
		return strings.Join(regions1, ",")
	}
	if len(regions2) > 0 {
		return strings.Join(regions2, ",")
	}
	if len(regions3) > 0 {
		return strings.Join(regions3, ",")
	}
	if len(regions4) > 0 {
		return strings.Join(regions4, ",")
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

type TransAnnoRow struct {
	Snv    string
	Gene   string
	GeneID string
	Region string
	Event  string
	Detail string
}

func (this TransAnnoRow) PK() string {
	return fmt.Sprintf("%s\t%s", this.Snv, this.Gene)
}

type TransAnnoScanner struct {
	pkg.IOScanner
	FieldNames []string
	DBname     string
}

func NewTransAnnoScanner(reader io.ReadCloser, dbname string) TransAnnoScanner {
	scanner := pkg.NewIOScanner(reader)
	scanner.Scan()
	return TransAnnoScanner{IOScanner: scanner, FieldNames: strings.Split(scanner.Text(), "\t"), DBname: dbname}
}

func (this TransAnnoScanner) Row() (TransAnnoRow, error) {
	row := make(map[string]string)
	fields := strings.Split(this.Text(), "\t")
	for i, fieldName := range this.FieldNames {
		row[fieldName] = fields[i]
	}
	return TransAnnoRow{
		Snv:    fmt.Sprintf("%s\t%s\t%s\t%s\t%s", row["Chr"], row["Start"], row["End"], row["Ref"], row["Alt"]),
		Gene:   row[this.DBname+".Gene"],
		GeneID: row[this.DBname+".GeneID"],
		Event:  row[this.DBname+".Event"],
		Region: row[this.DBname+".Region"],
		Detail: row[this.DBname+".Detail"],
	}, nil
}

func AggsTransAnno(transAnnoOutput, dbname, outfile string) error {
	reader, err := pkg.NewIOReader(transAnnoOutput)
	if err != nil {
		return err
	}
	defer reader.Close()
	scanner := NewTransAnnoScanner(reader, dbname)
	geneAnnos := make(map[string]GeneAnno)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return err
		}
		geneAnno, ok := geneAnnos[row.PK()]
		if !ok {
			geneAnno = GeneAnno{
				Snv:    row.Snv,
				Gene:   row.Gene,
				GeneID: row.GeneID,
			}
		}
		geneAnno.AddInfo(row.Event, row.Region, row.Detail)
		geneAnnos[row.PK()] = geneAnno
	}
	geneAnnoTexts := make([]string, len(geneAnnos))
	i := 0
	for _, snvGeneAnno := range geneAnnos {
		geneAnnoTexts[i] = snvGeneAnno.Text()
		i++
	}
	sort.Strings(geneAnnoTexts)
	writer, err := os.Create(outfile)
	if err != nil {
		return err
	}
	defer writer.Close()
	fmt.Fprintf(writer, "%s\n", strings.Join(scanner.FieldNames, "\t"))
	for _, text := range geneAnnoTexts {
		fmt.Fprintf(writer, "%s\n", text)
	}
	return nil
}

type AggsParam struct {
	AnnoGBOutput string `validate:"required,pathexists"`
	DBname       string `validate:"required"`
	Output       string `validate:"required"`
}

func (this AggsParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathexists", pkg.CheckPathExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return nil
}

func (this AggsParam) Run() error {
	return AggsTransAnno(this.AnnoGBOutput, this.DBname, this.Output)
}

func NewAggsCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "aggs",
		Short: "Aggregate Snv GeneBased Annotation",
		Run: func(cmd *cobra.Command, args []string) {
			var param AggsParam
			param.AnnoGBOutput, _ = cmd.Flags().GetString("gbanno")
			param.DBname, _ = cmd.Flags().GetString("dbname")
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
	cmd.Flags().StringP("gbanno", "g", "", "Annotate Result File of GeneBased")
	cmd.Flags().StringP("dbname", "n", "", "DB Name of GeneBased")
	cmd.Flags().StringP("output", "o", "", "Output Aggregated Annotation File")
	return cmd
}
