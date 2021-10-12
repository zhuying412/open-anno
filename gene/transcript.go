package gene

import "grandanno/db"

type Transcript struct {
	Chrom      string `json:"chrom"`
	Start      int    `json:"start"`
	End        int    `json:"end"`
	Transcript string `json:"transcript"`
	Sequence   string `json:"sequence"`
}

func (t Transcript) UpSteam() int {
	return t.Start - db.UpDownStreamLen
}

func (t Transcript) DownStream() int {
	return t.End + db.UpDownStreamLen
}
