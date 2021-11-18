package gene

type EntrezIdMap map[string]string

func (e EntrezIdMap) GetEntrezId(symbol string) string {
	if entrezId, ok := e[symbol]; ok {
		return entrezId
	}
	return ""
}
