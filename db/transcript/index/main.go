package index

import (
	"OpenAnno/db/chromosome"
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

func ReadRefgeneLine(refgeneLine string, upDownStreamLen int) Transcript {
	trans := Transcript{}
	fields := strings.Split(refgeneLine, "\t")
	trans.ID = fields[1]
	trans.Chrom = strings.Replace(strings.Split(fields[2], " ")[0], "chr", "", -1)
	if start, err := strconv.Atoi(fields[4]); err == nil {
		trans.Start = start + 1
	} else {
		log.Panic(err)
	}
	if end, err := strconv.Atoi(fields[5]); err == nil {
		trans.End = end
	} else {
		log.Panic(err)
	}
	trans.UpStream = trans.Start - upDownStreamLen
	trans.DownStream = trans.End + upDownStreamLen
	return trans
}

func ReadRefgeneFile(RefgeneFile string, upDownStreamLen int) Transcripts {
	fi, err := os.Open(RefgeneFile)
	defer func(fp *os.File) {
		err = fp.Close()
		if err != nil {
			log.Panic(err)
		}
	}(fi)
	if err != nil {
		log.Panic(err)
	}
	transcripts := make(Transcripts, 0)
	reader := bufio.NewReader(fi)
	for {
		line, err := reader.ReadBytes('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Panic(err)
			}
		}
		trans := ReadRefgeneLine(string(line), upDownStreamLen)
		if trans.Chrom == "M" || len(trans.Chrom) > 2 {
			continue
		}
		transcripts = append(transcripts, trans)
	}
	sort.Sort(transcripts)
	return transcripts
}

func InitTranscriptIndexes(indexStepLen int) TranscriptIndexes {
	indexes := make(TranscriptIndexes, 0)
	for _, chrom := range chromosome.ChromList {
		for i := 0; i < chrom.Length; i += indexStepLen {
			end := i + indexStepLen
			if end > chrom.Length {
				end = chrom.Length
			}
			index := TranscriptIndex{Chrom: chrom.Name, Start: i + 1, End: end}
			indexes = append(indexes, index)
		}
	}
	sort.Sort(indexes)
	return indexes
}

func ReadTranscriptIndexJSON(transcriptIndexFile string) TranscriptIndexes {
	fi, err := os.Open(transcriptIndexFile)
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
	indexes := make(TranscriptIndexes, 0)
	for {
		line, err := reader.ReadBytes('\n')
		if err == nil {
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
		var index TranscriptIndex
		err = json.Unmarshal(line, &index)
		if err != nil {
			log.Panic(err)
		}
		indexes = append(indexes, index)
	}
	sort.Sort(indexes)
	return indexes
}

func CreateTranscriptIndexJSON(indexes TranscriptIndexes, transcripts Transcripts, transcriptIndexFile string) {
	fo, err := os.Create(transcriptIndexFile)
	if err != nil {
		log.Panic(err)
	}
	defer func(fo *os.File) {
		err = fo.Close()
		if err != nil {
			log.Panic(err)
		}
	}(fo)
	for _, index := range indexes {
		index.SetTranscript(transcripts)
		if len(index.Transcripts) > 0 {
			contents, err := json.Marshal(index)
			if err != nil {
				log.Panic(err)
			}
			_, err = fo.Write(contents)
			if err != nil {
				log.Panic(err)
			}
		}
	}
}
