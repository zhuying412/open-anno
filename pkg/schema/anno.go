package schema

import (
	"fmt"
	"open-anno/pkg"
	"strings"
)

type SnvGeneAnno struct {
	Snv     string
	Gene    string
	GeneID  string
	Regions []string
	Events  []string
	Details []string
}

func (this *SnvGeneAnno) AddInfo(event, region, detail string) {
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

func (this SnvGeneAnno) Region() string {
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

func (this SnvGeneAnno) Event() string {
	if len(this.Events) > 0 {
		return strings.Join(this.Events, ",")
	}
	return "."
}

func (this SnvGeneAnno) Detail() string {
	if len(this.Details) > 0 {
		return strings.Join(this.Details, ",")
	}
	return "."
}

func (this SnvGeneAnno) Text() string {
	return fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s", this.Snv, this.Gene, this.GeneID, this.Event(), this.Region(), this.Detail())
}
