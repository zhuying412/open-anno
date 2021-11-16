package transcript

import (
	"OpenAnno/db/gene"
	"OpenAnno/seq"
	"bufio"
	"bytes"
	"encoding/json"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

func NewPosition(refgeneLineField []string) Position {
	position := Position{}
	var err error
	if position.ExonStart, err = strconv.Atoi(refgeneLineField[4]); err != nil {
		log.Panic(err)
	} else {
		position.ExonStart += 1
	}
	if position.ExonEnd, err = strconv.Atoi(refgeneLineField[5]); err != nil {
		log.Panic(err)
	}
	if position.CdsStart, err = strconv.Atoi(refgeneLineField[6]); err != nil {
		log.Panic(err)
	} else {
		position.CdsStart += 1
	}
	if position.CdsEnd, err = strconv.Atoi(refgeneLineField[7]); err != nil {
		log.Panic(err)
	}
	for _, strPos := range strings.Split(strings.Trim(refgeneLineField[9], ","), ",") {
		if pos, err := strconv.Atoi(strPos); err == nil {
			position.ExonStarts = append(position.ExonStarts, pos+1)
		}
	}
	for _, strPos := range strings.Split(strings.Trim(refgeneLineField[10], ","), ",") {
		if pos, err := strconv.Atoi(strPos); err == nil {
			position.ExonEnds = append(position.ExonEnds, pos)
		}
	}
	return position
}

func ReadRefgeneLine(refgeneLine string, upDownStreamLen int) Transcript {
	field := strings.Split(refgeneLine, "\t")
	trans := Transcript{
		Position:   NewPosition(field),
		Chrom:      strings.Replace(field[2], "chr", "", 1),
		Transcript: field[1],
		Strand:     field[3][0],
		Gene:       field[12],
		Tag:        field[13],
		EntrezId:   gene.SymbolToEntrezId.GetEntrezId(field[12]),
	}
	exonNum := len(trans.ExonStarts)
	for i := 0; i < exonNum; i++ {
		if i > 0 {
			region := Region{
				Start: trans.ExonEnds[i-1] + 1,
				End:   trans.ExonStarts[i] - 1,
				Type:  "intron",
			}
			trans.Regions = append(trans.Regions, region)
		}
		var exonOrder int
		if trans.Strand == '+' {
			exonOrder = i + 1
		} else {
			exonOrder = exonNum - 1
		}
		start, end := trans.ExonStarts[i], trans.ExonEnds[i]
		if trans.CdsStart > end || trans.CdsEnd < start ||
			trans.CdsStart <= start && end <= trans.CdsEnd {
			var typo string
			if trans.CdsStart > end {
				if trans.Strand == '+' {
					typo = "UTR5"
				} else {
					typo = "UTR3"
				}
			} else if trans.CdsEnd < start {
				if trans.Strand == '-' {
					typo = "UTR5"
				} else {
					typo = "UTR3"
				}
			} else {
				typo = "CDS"
			}
			trans.Regions = append(trans.Regions, Region{
				Start:     start,
				End:       end,
				Type:      typo,
				ExonOrder: exonOrder,
			})
		} else {
			UTRType1, UTRType2 := "", ""
			cdsStart, cdsEnd := start, end
			if start <= trans.CdsStart && trans.CdsStart <= end {
				if trans.Strand == '+' {
					UTRType1 = "UTR5"
				} else {
					UTRType1 = "UTR3"
				}
				cdsStart = trans.CdsStart
			}
			if start <= trans.CdsEnd && trans.CdsEnd <= end {
				if trans.Strand == '+' {
					UTRType2 = "UTR3"
				} else {
					UTRType2 = "UTR5"
				}
				cdsEnd = trans.CdsEnd
			}
			if UTRType1 != "" {
				trans.Regions = append(trans.Regions, Region{
					Start:     start,
					End:       trans.CdsStart - 1,
					Type:      UTRType1,
					ExonOrder: exonOrder,
				})
			}
			trans.Regions = append(trans.Regions, Region{
				Start:     cdsStart,
				End:       cdsEnd,
				Type:      "CDS",
				ExonOrder: exonOrder,
			})
			if UTRType2 != "" {
				trans.Regions = append(trans.Regions, Region{
					Start:     trans.CdsEnd + 1,
					End:       end,
					Type:      UTRType2,
					ExonOrder: exonOrder,
				})
			}
		}
	}
	sort.Sort(trans.Regions)
	stream1 := Region{
		Start: trans.ExonStart - upDownStreamLen,
		End:   trans.ExonStart - 1,
	}
	stream2 := Region{
		Start: trans.ExonEnd + 1,
		End:   trans.ExonEnd + upDownStreamLen,
	}
	if trans.Strand == '+' {
		stream1.Type = "upstream"
		stream2.Type = "downstream"
	} else {
		stream1.Type = "upstream"
		stream2.Type = "downstream"
	}
	trans.Streams = Regions{stream1, stream2}
	return trans
}

func ReadRefgeneFile(refgeneFile string, upDownStreamLen int) Transcripts {
	fi, err := os.Open(refgeneFile)
	if err != nil {
		log.Panic(err)
	}
	defer func(fp *os.File) {
		err := fp.Close()
		if err != nil {
			log.Panic(err)
		}
	}(fi)
	transcripts := make(Transcripts, 0)
	reader := bufio.NewReader(fi)
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Panic(err)
			}
		}
		line = strings.TrimSpace(line)
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		trans := ReadRefgeneLine(line, upDownStreamLen)
		if trans.Chrom == "M" || len(trans.Chrom) > 2 {
			continue
		}
		transcripts = append(transcripts, trans)
	}
	sort.Sort(transcripts)
	return transcripts
}

func ReadTranscriptJSON(transcriptFile string) TranscriptMap {
	fi, err := os.Open(transcriptFile)
	if err != nil {
		log.Panic(err)
	}
	defer func(fp *os.File) {
		err := fp.Close()
		if err != nil {
			log.Panic(err)
		}
	}(fi)
	reader := bufio.NewReader(fi)
	transMap := make(TranscriptMap)
	for {
		line, err := reader.ReadBytes('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Panic(err)
			}
		}
		line = bytes.TrimSpace(line)
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		var trans Transcript
		err = json.Unmarshal(line, &trans)
		if err != nil {
			log.Panic(err)
		}
		transMap[trans.SN()] = trans
	}
	return transMap
}

func CreateTranscriptJSON(transcripts []Transcript, chromSeq seq.Sequence, transcriptFile string) {
	fo, err := os.Create(transcriptFile)
	if err != nil {
		log.Panic(err)
	}
	defer func(fo *os.File) {
		err = fo.Close()
		if err != nil {
			log.Panic(err)
		}
	}(fo)
	for _, trans := range transcripts {
		sequence := chromSeq.SubSeq(trans.ExonStart-1, trans.ExonEnd-trans.ExonStart+1)
		trans.SetSequence(sequence)
		contents, err := json.Marshal(trans)
		if err != nil {
			log.Panic(err)
		}
		_, err = fo.Write(contents)
		if err != nil {
			log.Panic(err)
		}
		_, err = fo.Write([]byte{'\n'})
		if err != nil {
			log.Panic(err)
		}
	}
}
