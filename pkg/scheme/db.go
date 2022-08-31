package scheme

type DBVar struct {
	Chrom string `json:"chrom"`
	Start int    `json:"start"`
	End   int    `json:"end"`
	Ref   string `json:"ref"`
	Alt   string `json:"alt"`
	Text  string `json:"text"`
}

func (this DBVar) Compare(variant Variant) string {
	dbvar := Variant{
		Chrom: this.Chrom,
		Start: this.Start,
		End:   this.End,
		Ref:   this.Ref,
		Alt:   this.Alt,
	}
	return dbvar.Compare(variant)
}

type DBVarIdx struct {
	Bin   string `json:"bin"`
	Start int64  `json:"start"`
	End   int64  `json:"end"`
}

type DBReg struct {
	Chrom string `json:"chrom"`
	Start int    `json:"start"`
	End   int    `json:"end"`
	Info  string `json:"info"`
}
