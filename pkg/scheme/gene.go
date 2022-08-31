package scheme

type GeneInfo struct {
	EntrezId string `json:"entrez_id"`
	Symbol   string `json:"symbol"`
	Chrom    string `json:"chrom"`
}

type GeneInfoMap map[string]map[string]GeneInfo
