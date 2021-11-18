package cnv

import (
	"OpenAnno/pkg/transcript"
	"OpenAnno/pkg/transcript/index"
	"OpenAnno/pkg/variant"
)

func NewGeneAnnoInIntergeic() GeneAnno {
	return GeneAnno{{Region: "intergenic"}}
}

func NewGeneAnnoInUpDownStream(cnv variant.Cnv, transcripts []transcript.Transcript) GeneAnno {
	annos := make(GeneAnno, 0)
	for _, trans := range transcripts {
		for _, region := range trans.Streams {
			if cnv.Start <= region.End && cnv.End >= region.Start {
				var _anno GeneAnnoItem
				_anno.SetGene(trans)
				_anno.SetRegion(region.Type)
				annos = append(annos, _anno)
				break
			}
		}
	}
	return annos
}

func NewGeneAnnoInGene(cnv variant.Cnv, transcripts []transcript.Transcript) GeneAnno {
	var cmplAnnos, incmplAnnos, unkAnnos, annos GeneAnno
	for _, trans := range transcripts {
		if cnv.End >= trans.ExonStart && cnv.Start <= trans.ExonEnd {
			var _anno GeneAnnoItem
			_anno.SetGene(trans)
			if trans.IsUnk() {
				_anno.SetRegion("ncRNA")
				unkAnnos = append(unkAnnos, _anno)
			} else {
				for _, region := range trans.Regions {
					if cnv.Start <= region.End && cnv.End >= region.Start {
						if region.IsCDS() {
							_anno.SetExon(region.ExonOrder)
						}
					}
				}
				if trans.IsCmpl() {
					_anno.SetRegion("Gene")
					cmplAnnos = append(cmplAnnos, _anno)
				} else {
					_anno.SetRegion("incmplCDS")
					incmplAnnos = append(incmplAnnos, _anno)
				}
			}
		}
	}
	annos = append(annos, cmplAnnos...)
	if !annos.HasCodingorSplicingRegion() {
		annos = append(annos, incmplAnnos...)
	}
	if len(annos) == 0 {
		annos = append(annos, unkAnnos...)
	}
	return annos
}

func RunAnnotate(cnvs variant.Cnvs, transMap transcript.TranscriptMap, transIndexes index.TranscriptIndexes) map[string]GeneAnno {
	annoMap := make(map[string]GeneAnno)
	for i := 0; i < len(cnvs); i++ {
		for j := 0; j < len(transIndexes); j++ {
			refStart, refEnd := transIndexes[j].Start, transIndexes[j].End
			if cnvs[i].End < refStart {
				break
			} else if cnvs[i].Start > refEnd {
				continue
			} else {
				var annos GeneAnno
				if cnvs[i].End < refStart {
					annos = NewGeneAnnoInIntergeic()
				} else {
					transcripts := transMap.FindMany(transIndexes[j].Transcripts)
					annos = NewGeneAnnoInGene(cnvs[i], transcripts)
					if len(annos) == 0 {
						annos = NewGeneAnnoInUpDownStream(cnvs[i], transcripts)
					}
					if len(annos) == 0 {
						annos = NewGeneAnnoInIntergeic()
					}
				}
				annoMap[cnvs[i].SN()] = annos
			}
		}
	}
	return annoMap
}
