package anno

type AnnoType string

const (
	AnnoType_GENE AnnoType = "gene_based"
	AnnoType_DB   AnnoType = "database_based"
)

type IAnno interface {
	AnnoType() AnnoType
}
