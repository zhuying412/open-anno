package anno

type AnnoType string

const (
	AnnoType_GENE   AnnoType = "gene_based"
	AnnoType_FILTER AnnoType = "filter_based"
	AnnoType_REGION AnnoType = "region_based"
)

type IAnno interface {
	AnnoType() AnnoType
}
