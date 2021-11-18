package del

import (
	"OpenAnno/pkg/transcript"
	"sort"
	"strings"
)

func (a *GeneAnnoItem) AnnoInMulti(regions transcript.Regions) {
	regionTypes := make([]string, 0)
	for _, region := range regions {
		if sort.SearchStrings(regionTypes, string(region.Type)) != -1 {
			regionTypes = append(regionTypes, string(region.Type))
		}
	}
	sort.Strings(regionTypes)
	a.SetRegion(strings.Join(regionTypes, "_"))
}
