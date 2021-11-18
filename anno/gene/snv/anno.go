package snv

import (
	"OpenAnno/anno"
	"OpenAnno/anno/gene/snv/del"
	"OpenAnno/anno/gene/snv/ins"
	"OpenAnno/anno/gene/snv/snp"
	"OpenAnno/pkg/transcript"
	"OpenAnno/pkg/transcript/index"
	"OpenAnno/pkg/variant"
)

type IGeneAnnoItem interface {
	InCodingOrSplicingRegion() bool
	SetExon(exon int)
	SetGene(trans transcript.Transcript)
	SetRegion(region string)
	SetEvent(event string)
	SetNAChange(change string)
	SetAAChange(change string)
	AnnoInGene(snv variant.Snv, trans transcript.Transcript)
}

type GeneAnno []IGeneAnnoItem

func (a GeneAnno) AnnoType() anno.AnnoType {
	return anno.AnnoType_GENE
}

func (a GeneAnno) HasCodingorSplicingRegion() bool {
	for _, _anno := range a {
		if _anno.InCodingOrSplicingRegion() {
			return true
		}
	}
	return false
}

func NewGeneAnnoItem(snvType variant.SnvType) IGeneAnnoItem {
	if snvType == variant.SnvType_DEL {
		return &del.GeneAnnoItem{}
	} else if snvType == variant.SnvType_INS {
		return &ins.GeneAnnoItem{}
	} else {
		return &snp.GeneAnnoItem{}
	}
}

func NewGeneAnnoInIntergeic(snvType variant.SnvType) GeneAnno {
	annoItem := NewGeneAnnoItem(snvType)
	annoItem.SetRegion("intergenic")
	return GeneAnno{annoItem}
}

func NewGeneAnnosInUpDownStream(snv variant.Snv, transcripts transcript.Transcripts) GeneAnno {
	_anno := make(GeneAnno, 0)
	for _, trans := range transcripts {
		for _, region := range trans.Streams {
			if snv.Start <= region.End && snv.End >= region.Start {
				annoItem := NewGeneAnnoItem(snv.Type())
				annoItem.SetGene(trans)
				annoItem.SetRegion(string(region.Type))
				_anno = append(_anno, annoItem)
				break
			}
		}
	}
	return _anno
}

func NewGeneAnnosInGene(snv variant.Snv, transcripts transcript.Transcripts) GeneAnno {
	var _anno, cmplAnno, incmplAnno, unkAnno GeneAnno
	for _, trans := range transcripts {
		if snv.End >= trans.ExonStart && snv.Start <= trans.ExonEnd {
			annoItem := NewGeneAnnoItem(snv.Type())
			if trans.IsUnk() {
				annoItem.SetGene(trans)
				annoItem.SetRegion("ncRNA")
				unkAnno = append(unkAnno, annoItem)
			} else {
				annoItem.AnnoInGene(snv, trans)
				if trans.IsCmpl() {
					cmplAnno = append(cmplAnno, annoItem)
				} else {
					annoItem.SetEvent("incmplCDS")
					incmplAnno = append(incmplAnno, annoItem)
				}
			}
		}
	}
	_anno = append(_anno, cmplAnno...)
	if !_anno.HasCodingorSplicingRegion() {
		_anno = append(_anno, incmplAnno...)
	}
	if len(_anno) == 0 {
		_anno = append(_anno, unkAnno...)
	}
	return _anno
}

func RunAnnotate(snvs variant.Snvs, transMap transcript.TranscriptMap, transIndexes index.TranscriptIndexes) map[string]GeneAnno {
	snp.Init()
	annoMap := make(map[string]GeneAnno)
	for i, j := 0, 0; i < len(snvs) && j < len(transIndexes); {
		if snvs[i].Start > transIndexes[j].End {
			j++
		} else {
			var annos GeneAnno
			if snvs[i].End < transIndexes[j].Start {
				annos = NewGeneAnnoInIntergeic(snvs[i].Type())
			} else {
				transcripts := transMap.FindMany(transIndexes[j].Transcripts)
				annos = NewGeneAnnosInGene(snvs[i], transcripts)
				if len(annos) == 0 {
					annos = NewGeneAnnosInUpDownStream(snvs[i], transcripts)
				}
				if len(annos) == 0 {
					annos = NewGeneAnnoInIntergeic(snvs[i].Type())
				}
			}
			annoMap[snvs[i].SN()] = annos
			i++
		}
	}
	return annoMap
}