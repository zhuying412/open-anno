package transcript

import (
	"OpenAnno/pkg/seq"
	"fmt"
)

type TransTag string

const (
	TransTag_CMPL   TransTag = "cmpl"
	TransTag_UNK    TransTag = "unk"
	TransTag_INCMPL TransTag = "incmpl"
)

type Position struct {
	ExonStart  int   `json:"exon_start"`
	ExonEnd    int   `json:"exon_end"`
	CdsStart   int   `json:"cds_start"`
	CdsEnd     int   `json:"cds_end"`
	ExonStarts []int `json:"exon_starts"`
	ExonEnds   []int `json:"exon_ends"`
}

type Transcript struct {
	Chrom      string `json:"chrom"`
	Strand     byte   `json:"strand"`
	Gene       string `json:"gene"`
	EntrezId   string `json:"entrez_id"`
	Transcript string `json:"transcript"`
	Position   `json:"position"`
	Regions    Regions      `json:"regions"`
	Streams    Regions      `json:"streams"`
	Tag        TransTag     `json:"tag"`
	Mrna       seq.Sequence `json:"mrna"`
	Cdna       seq.Sequence `json:"cdna"`
	Protein    seq.Sequence `json:"protein"`
}

func (t Transcript) SN() string {
	return fmt.Sprintf("%s|%s:%d:%d", t.Transcript, t.Chrom, t.ExonStart, t.ExonEnd)
}

func (t Transcript) IsCmpl() bool {
	return t.Tag == TransTag_CMPL
}

func (t Transcript) IsUnk() bool {
	return t.Tag == TransTag_UNK
}

func (t Transcript) InChromMT() bool {
	return t.Chrom == "MT"
}

func (t *Transcript) SetSequence(sequence seq.Sequence) {
	if !sequence.IsEmpty() {
		t.Mrna = sequence
		if t.IsUnk() {
			for _, region := range t.Regions {
				if region.IsCDS() {
					t.Cdna.Join(t.Mrna.SubSeq(region.Start-t.ExonStart, region.End-region.Start+1))
				}
			}
		}
		if t.Strand == '-' {
			t.Mrna.ReverseComplementing()
			t.Cdna.ReverseComplementing()
		}
		if !t.Cdna.IsEmpty() {
			t.Protein = t.Cdna.Translate(t.InChromMT())
			if t.Protein.IsProteinCmpl() {
				t.Tag = TransTag_CMPL
			} else {
				t.Tag = TransTag_INCMPL
				t.Protein.Clear()
			}
		} else {
			t.Tag = TransTag_UNK
		}
	}
}

type Transcripts []Transcript

func (Transcripts Transcripts) Len() int {
	return len(Transcripts)
}

func (t Transcripts) Less(i, j int) bool {
	return t[i].ExonStart < t[j].ExonStart || t[i].ExonEnd < t[j].ExonEnd
}

func (t Transcripts) Swap(i, j int) {
	t[i], t[j] = t[j], t[i]
}

func (t Transcripts) FilterByChrom(chrom string) Transcripts {
	transcripts := make(Transcripts, 0)
	for _, trans := range t {
		if trans.Chrom == chrom {
			transcripts = append(transcripts, trans)
		}
	}
	return transcripts
}

type TranscriptMap map[string]Transcript

func (t TranscriptMap) FilterByChrom(chrom string) TranscriptMap {
	transMap := make(TranscriptMap)
	for sn, trans := range t {
		if trans.Chrom == chrom {
			transMap[sn] = trans
		}
	}
	return transMap
}

func (t TranscriptMap) FindMany(sns []string) Transcripts {
	transcripts := make(Transcripts, 0)
	for _, sn := range sns {
		if trans, ok := t[sn]; ok {
			transcripts = append(transcripts, trans)
		}
	}
	return transcripts
}
