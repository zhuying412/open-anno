package gene

import (
	"bufio"
	"fmt"
	"grandanno/db"
	"grandanno/seq"
	"io"
	"log"
	"os"
	"sort"
	"strings"
)

type Refgene struct {
	Chrom      string       `json:"chrom"`
	Strand     byte         `json:"strand"`
	Gene       string       `json:"gene"`
	EntrezId   string       `json:"entrez_id"`
	Transcript string       `json:"transcript"`
	Position   Position     `json:"position"`
	Regions    Regions      `json:"regions"`
	Streams    Regions      `json:"streams"`
	Tag        string       `json:"tag"`
	Mrna       seq.Sequence `json:"mrna"`
	Cdna       seq.Sequence `json:"cdna"`
	Protein    seq.Sequence `json:"protein"`
}

func (r Refgene) SN() string {
	return fmt.Sprintf("%s|%s:%d:%d", r.Transcript, r.Chrom, r.Position.ExonStart, r.Position.ExonEnd)
}

func (r Refgene) IsCmpl() bool {
	return r.Tag == "cmpl"
}

func (r Refgene) IsUnk() bool {
	return r.Tag == "unk"
}

func (r *Refgene) SetSequence(sequence seq.Sequence) {
	if !sequence.IsEmpty() {
		r.Mrna = sequence
		if r.Tag != "unk" {
			for _, region := range r.Regions {
				if region.Type == "CDS" {
					r.Cdna.Join(r.Mrna.SubSeq(region.Start-r.Position.ExonStart, region.End-region.Start+1))
				}
			}
			if r.Strand == '-' {
				r.Cdna.Reverse()
			}
		}
		if !r.Cdna.IsEmpty() {
			r.Protein = r.Cdna.Translate(r.Chrom == "MT")
			if r.Protein.IsProteinCmpl() {
				r.Tag = "cmpl"
			} else {
				r.Tag = "incmpl"
				r.Protein.Clear()
			}
		} else {
			r.Tag = "unk"
		}
	}
}

func NewRefgene(refgeneLine string) Refgene {
	field := strings.Split(refgeneLine, "\t")
	refgene := Refgene{
		Position:   NewPosition(field),
		Chrom:      strings.Replace(field[2], "chr", "", 1),
		Transcript: field[1],
		Strand:     field[3][0],
		Gene:       field[12],
		Tag:        field[13],
		EntrezId:   db.SymbolToEntrezId.GetEntrezId(field[12]),
	}
	exonNum := len(refgene.Position.ExonStarts)
	for i := 0; i < exonNum; i++ {
		if i > 0 {
			region := Region{
				Start: refgene.Position.ExonEnds[i-1] + 1,
				End:   refgene.Position.ExonStarts[i] - 1,
				Type:  "intron",
			}
			refgene.Regions = append(refgene.Regions, region)
		}
		var exonOrder int
		if refgene.Strand == '+' {
			exonOrder = i + 1
		} else {
			exonOrder = exonNum - 1
		}
		start, end := refgene.Position.ExonStarts[i], refgene.Position.ExonEnds[i]
		if refgene.Position.CdsStart > end || refgene.Position.CdsEnd < start ||
			refgene.Position.CdsStart <= start && end <= refgene.Position.CdsEnd {
			var typo string
			if refgene.Position.CdsStart > end {
				if refgene.Strand == '+' {
					typo = "UTR5"
				} else {
					typo = "UTR3"
				}
			} else if refgene.Position.CdsEnd < start {
				if refgene.Strand == '-' {
					typo = "UTR5"
				} else {
					typo = "UTR3"
				}
			} else {
				typo = "CDS"
			}
			refgene.Regions = append(refgene.Regions, Region{
				Start:     start,
				End:       end,
				Type:      typo,
				ExonOrder: exonOrder,
			})
		} else {
			UTRType1, UTRType2 := "", ""
			cdsStart, cdsEnd := start, end
			if start < refgene.Position.CdsStart && refgene.Position.CdsStart < end {
				if refgene.Strand == '+' {
					UTRType1 = "UTR5"
				} else {
					UTRType1 = "UTR3"
				}
				cdsStart = refgene.Position.CdsStart
			}
			if start < refgene.Position.CdsEnd && refgene.Position.CdsEnd < end {
				if refgene.Strand == '+' {
					UTRType2 = "UTR3"
				} else {
					UTRType2 = "UTR5"
				}
				cdsEnd = refgene.Position.CdsEnd
			}
			if UTRType1 != "" {
				refgene.Regions = append(refgene.Regions, Region{
					Start:     start,
					End:       refgene.Position.CdsStart - 1,
					Type:      UTRType1,
					ExonOrder: exonOrder,
				})
			}
			refgene.Regions = append(refgene.Regions, Region{
				Start:     cdsStart,
				End:       cdsEnd,
				Type:      "CDS",
				ExonOrder: exonOrder,
			})
			if UTRType2 != "" {
				refgene.Regions = append(refgene.Regions, Region{
					Start:     refgene.Position.CdsEnd + 1,
					End:       end,
					Type:      UTRType2,
					ExonOrder: exonOrder,
				})
			}
		}
	}
	sort.Sort(refgene.Regions)
	stream1 := Region{
		Start: refgene.Position.ExonStart - db.UpDownStreamLen,
		End:   refgene.Position.ExonStart - 1,
	}
	stream2 := Region{
		Start: refgene.Position.ExonEnd + 1,
		End:   refgene.Position.ExonEnd + db.UpDownStreamLen,
	}
	if refgene.Strand == '+' {
		stream1.Type = "upstream"
		stream2.Type = "downstream"
	} else {
		stream1.Type = "upstream"
		stream2.Type = "downstream"
	}
	refgene.Streams = Regions{stream1, stream2}
	return refgene
}

type RefgeneMap map[string]Refgene

func ReadRefgeneFile(refgeneFile string) RefgeneMap {
	refgeneMap := make(RefgeneMap)
	if fp, err := os.Open(refgeneFile); err == nil {
		defer func(fp *os.File) {
			err := fp.Close()
			if err != nil {
				log.Panic(err)
			}
		}(fp)
		reader := bufio.NewReader(fp)
		for {
			if line, err := reader.ReadString('\n'); err == nil {
				line = strings.TrimSpace(line)
				if len(line) == 0 || line[0] == '#' {
					continue
				}
				refgene := NewRefgene(line)
				if refgene.Chrom == "M" || len(refgene.Chrom) > 2 {
					continue
				}
				refgeneMap[refgene.SN()] = refgene
			} else {
				if err == io.EOF {
					break
				} else {
					log.Panic(err)
				}
			}
		}
	} else {
		log.Panic(err)
	}
	return refgeneMap
}

func (r *RefgeneMap) SetSequence(mrna seq.Fasta) {
	for sn, refgene := range *r {
		if sequence, ok := mrna[sn]; ok {
			refgene.SetSequence(sequence)
			(*r)[sn] = refgene
		}
	}
}

func (r RefgeneMap) FilterByChrom(chrom string) RefgeneMap {
	refgeneMap := make(RefgeneMap)
	for sn, refgene := range r {
		if refgene.Chrom == chrom {
			refgeneMap[sn] = refgene
		}
	}
	return refgeneMap
}

func (r RefgeneMap) FindMany(sns []string) (refgenes []Refgene) {
	refgenes = make([]Refgene, 0)
	for _, sn := range sns {
		if refgene, ok := r[sn]; ok {
			refgenes = append(refgenes, refgene)
		}
	}
	return refgenes
}
