package transcript

import (
	"OpenAnno/db"
	"OpenAnno/pkg/seq"
	"OpenAnno/pkg/transcript"
	"OpenAnno/pkg/utils"
	"log"
	"os"
	"path"
	"sort"
	"strings"
)

func newPosition(refgeneLineField []string) transcript.Position {
	position := transcript.Position{
		ExonStart: utils.StrToInt(refgeneLineField[4]) + 1,
		ExonEnd:   utils.StrToInt(refgeneLineField[5]) + 1,
		CdsStart:  utils.StrToInt(refgeneLineField[6]),
		CdsEnd:    utils.StrToInt(refgeneLineField[7]),
	}
	strExonStarts := strings.Split(strings.Trim(refgeneLineField[9], ","), ",")
	strExonEnds := strings.Split(strings.Trim(refgeneLineField[10], ","), ",")
	for i := 0; i < len(strExonStarts); i++ {
		position.ExonStarts = append(position.ExonStarts, utils.StrToInt(strExonStarts[i])+1)
		position.ExonEnds = append(position.ExonEnds, utils.StrToInt(strExonEnds[i]))
	}
	return position
}

func readRefgeneLine(refgeneLine string, upDownStreamLen int) transcript.Transcript {
	field := strings.Split(refgeneLine, "\t")
	trans := transcript.Transcript{
		Position:   newPosition(field),
		Chrom:      strings.Replace(field[2], "chr", "", 1),
		Transcript: field[1],
		Strand:     field[3][0],
		Gene:       field[12],
		Tag:        transcript.TransTag(field[13]),
		EntrezId:   db.SymbolToEntrezId.GetEntrezId(field[12]),
	}
	exonNum := len(trans.ExonStarts)
	for i := 0; i < exonNum; i++ {
		if i > 0 {
			region := transcript.Region{
				Start: trans.ExonEnds[i-1] + 1,
				End:   trans.ExonStarts[i] - 1,
				Type:  transcript.RegionType_INTRON,
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
			var typo transcript.RegionType
			if trans.CdsStart > end {
				if trans.Strand == '+' {
					typo = transcript.RegionType_UTR5
				} else {
					typo = transcript.RegionType_UTR3
				}
			} else if trans.CdsEnd < start {
				if trans.Strand == '-' {
					typo = transcript.RegionType_UTR5
				} else {
					typo = transcript.RegionType_UTR3
				}
			} else {
				typo = transcript.RegionType_CDS
			}
			trans.Regions = append(trans.Regions, transcript.Region{
				Start:     start,
				End:       end,
				Type:      typo,
				ExonOrder: exonOrder,
			})
		} else {
			var UTRType1, UTRType2 transcript.RegionType
			cdsStart, cdsEnd := start, end
			if start <= trans.CdsStart && trans.CdsStart <= end {
				if trans.Strand == '+' {
					UTRType1 = transcript.RegionType_UTR5
				} else {
					UTRType1 = transcript.RegionType_UTR3
				}
				cdsStart = trans.CdsStart
			}
			if start <= trans.CdsEnd && trans.CdsEnd <= end {
				if trans.Strand == '+' {
					UTRType2 = transcript.RegionType_UTR3
				} else {
					UTRType2 = transcript.RegionType_UTR5
				}
				cdsEnd = trans.CdsEnd
			}
			if UTRType1 != "" {
				trans.Regions = append(trans.Regions, transcript.Region{
					Start:     start,
					End:       trans.CdsStart - 1,
					Type:      UTRType1,
					ExonOrder: exonOrder,
				})
			}
			trans.Regions = append(trans.Regions, transcript.Region{
				Start:     cdsStart,
				End:       cdsEnd,
				Type:      "CDS",
				ExonOrder: exonOrder,
			})
			if UTRType2 != "" {
				trans.Regions = append(trans.Regions, transcript.Region{
					Start:     trans.CdsEnd + 1,
					End:       end,
					Type:      UTRType2,
					ExonOrder: exonOrder,
				})
			}
		}
	}
	sort.Sort(trans.Regions)
	stream1 := transcript.Region{
		Start: trans.ExonStart - upDownStreamLen,
		End:   trans.ExonStart - 1,
	}
	stream2 := transcript.Region{
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
	trans.Streams = transcript.Regions{stream1, stream2}
	return trans
}

func readRefgeneFile(refgeneFile string, upDownStreamLen int) transcript.Transcripts {
	fi, reader := utils.OpenFile(refgeneFile)
	defer utils.CloseFile(fi)
	transcripts := make(transcript.Transcripts, 0)
	for {
		line, isEof := utils.ReadLine(reader, '#')
		if isEof {
			break
		}
		trans := readRefgeneLine(line, upDownStreamLen)
		if trans.Chrom == "M" || len(trans.Chrom) > 2 {
			continue
		}
		transcripts = append(transcripts, trans)
	}
	sort.Sort(transcripts)
	return transcripts
}

func ReadTranscriptJSON(transcriptFile string) transcript.TranscriptMap {
	fi, reader := utils.OpenFile(transcriptFile)
	defer utils.CloseFile(fi)
	transMap := make(transcript.TranscriptMap)
	for {
		line, isEof := utils.ReadLine(reader, '#')
		if isEof {
			break
		}
		var trans transcript.Transcript
		utils.FromJSON(line, &trans)
		transMap[trans.SN()] = trans
	}
	return transMap
}

func CreateTranscriptJSON(refgeneFile string, referenceFile string, transcriptDir string, upDownStreamLen int) {
	utils.CreateDir(transcriptDir)
	log.Printf("read %s", referenceFile)
	reference := seq.ReadFastaFile(referenceFile)
	log.Printf("read %s", refgeneFile)
	transcripts := readRefgeneFile(refgeneFile, upDownStreamLen)
	for _, chrom := range db.ChromList {
		outfile := path.Join(transcriptDir, "chr"+chrom.Name+".json")
		_, err := os.Stat(outfile)
		if os.IsExist(err) {
			continue
		}
		log.Printf("write %s", outfile)
		subTranscripts := transcripts.FilterByChrom(chrom.Name)
		chromSeq := reference[chrom.Name]
		if subTranscripts.Len() > 0 {
			fo := utils.CreateFile(outfile)
			for _, trans := range subTranscripts {
				sequence := chromSeq.SubSeq(trans.ExonStart-1, trans.ExonEnd-trans.ExonStart+1)
				trans.SetSequence(sequence)
				contents := utils.ToJSON(trans)
				utils.WriteLine(fo, contents+"\n")
			}
			utils.CloseFile(fo)
		}
	}
}
