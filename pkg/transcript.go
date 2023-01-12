package pkg

import (
	"github.com/brentp/faidx"
)

// Transcript 转录本，继承自GenePred，加入GeneID和Regions信息
type Transcript struct {
	GenePred
	GeneID  string  `json:"gene_id"`
	Regions Regions `json:"regions"`
}

// HasUTR3 存在UTR3区域
func (this Transcript) HasUTR3() bool {
	if this.Strand == "+" {
		return this.CdsEnd < this.TxEnd
	}
	return this.CdsStart > this.TxStart
}

// HasUTR5 存在UTR5区域
func (this Transcript) HasUTR5() bool {
	if this.Strand == "+" {
		return this.CdsStart > this.TxStart
	}
	return this.CdsEnd < this.TxEnd
}

// CdsCount Regions中CDS元件数量
func (this Transcript) CdsCount() int {
	return this.Regions.CdsCount()
}

// CLen Regions中CDS区域总长度
func (this Transcript) CLen() int {
	return this.Regions.CLen()
}

// CDNA Regions中CDS的Sequence拼接结果即CDNA
func (this Transcript) CDNA() string {
	return this.Regions.CDNA()
}

// CDNA Regions中所有Sequence拼接结果即DNA
func (this Transcript) DNA() string {
	return this.Regions.DNA()
}

// SetGeneID 根据geneSymbolToID的Map信息设置转录本的GeneID
func (this *Transcript) SetGeneID(geneSymbolToID map[string]map[string]string) {
	if entrezId, ok := geneSymbolToID[this.Chrom][this.Gene]; ok {
		this.GeneID = entrezId
	}
}

// SetRegions 设置初始化Regions信息
func (this *Transcript) SetRegions() error {
	var err error
	this.Regions, err = NewRegions(*this)
	if err != nil {
		return err
	}
	return nil
}

// SetRegions 设置初始化Regions信息, 并设置每个Regions的Sequence
func (this *Transcript) SetRegionsWithSeq(genome *faidx.Faidx) error {
	var err error
	this.Regions, err = NewRegionsWithSeq(*this, genome)
	if err != nil {
		return err
	}
	return nil
}

// Transcripts 以PK为Key的Transcript的Map数据结构
type Transcripts map[string]Transcript

// NewTranscripts 创建指定染色体的Transcripts信息
func NewTranscripts(gpes GenePreds, chrom string, geneSymbolToID map[string]map[string]string) (Transcripts, error) {
	transcripts := make(Transcripts)
	for pk, gpe := range gpes {
		if gpe.Chrom == chrom {
			trans := Transcript{GenePred: gpe}
			trans.SetGeneID(geneSymbolToID)
			err := trans.SetRegions()
			if err != nil {
				return transcripts, err
			}
			transcripts[pk] = trans
		}
	}
	return transcripts, nil
}

// NewTranscriptsWithSeq 获取指定染色体的Transcripts信息，并设置Sequence
func NewTranscriptsWithSeq(gpes GenePreds, chrom string, geneSymbolToID map[string]map[string]string, genome *faidx.Faidx) (Transcripts, error) {
	transcripts := make(Transcripts)
	for sn, gpe := range gpes {
		if gpe.Chrom == chrom {
			trans := Transcript{GenePred: gpe}
			trans.SetGeneID(geneSymbolToID)
			err := trans.SetRegionsWithSeq(genome)
			if err != nil {
				return transcripts, err
			}
			transcripts[sn] = trans
		}
	}
	return transcripts, nil
}
