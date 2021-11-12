package region

import (
	"OpenAnno/anno"
	"OpenAnno/variant"
)

type RegionBasedAnno map[string]string

func (f RegionBasedAnno) AnnoType() anno.AnnoType {
	return anno.AnnoType_REGION
}

func ReadRegionFields(headers []string, fields []string) (variant.Variant, RegionBasedAnno) {
	_anno := make(RegionBasedAnno)
	commomVariant := variant.NewVariant(fields[0], fields[1], fields[2], fields[3], fields[4])
	for i := 5; i < len(fields); i++ {
		_anno[headers[i]] = fields[i]
	}
	return commomVariant, _anno
}
