package transcript

import (
	"OpenAnno/db"
	"OpenAnno/pkg/transcript/index"
	"OpenAnno/pkg/utils"
	"log"
	"os"
	"path"
	"sort"
	"strings"
)

func readRefgeneLineForIndex(refgeneLine string, upDownStreamLen int) index.Transcript {
	fields := strings.Split(refgeneLine, "\t")
	trans := index.Transcript{
		ID:    fields[1],
		Chrom: strings.Replace(strings.Split(fields[2], " ")[0], "chr", "", -1),
		Start: utils.StrToInt(fields[4]),
		End:   utils.StrToInt(fields[5]),
	}
	trans.UpStream = trans.Start - upDownStreamLen
	trans.DownStream = trans.End + upDownStreamLen
	return trans
}

func readRefgeneFileForIndex(refgeneFile string, upDownStreamLen int) index.Transcripts {
	fi, reader := utils.OpenFile(refgeneFile)
	defer utils.CloseFile(fi)
	transcripts := make(index.Transcripts, '#')
	for {
		line, isEof := utils.ReadLine(reader, 0)
		if isEof {
			break
		}
		trans := readRefgeneLineForIndex(line, upDownStreamLen)
		if trans.Chrom == "M" || len(trans.Chrom) > 2 {
			continue
		}
		transcripts = append(transcripts, trans)
	}
	sort.Sort(transcripts)
	return transcripts
}

func ReadTranscriptIndexJSON(transcriptIndexFile string) index.TranscriptIndexes {
	fi, reader := utils.OpenFile(transcriptIndexFile)
	defer utils.CloseFile(fi)
	indexes := make(index.TranscriptIndexes, 0)
	for {
		line, isEof := utils.ReadLine(reader, '#')
		if isEof {
			break
		}
		var _index index.TranscriptIndex
		utils.FromJSON(line, &_index)
		indexes = append(indexes, _index)
	}
	sort.Sort(indexes)
	return indexes
}

func CreateTranscriptIndexJSON(refgeneFile string, transcriptDir string, upDownStreamLen int, indexStepLen int) {
	utils.CreateDir(transcriptDir)
	log.Printf("read %s", refgeneFile)
	transcripts := readRefgeneFileForIndex(refgeneFile, upDownStreamLen)
	log.Print("init transcript indexes")
	indexes := index.InitTranscriptIndexes(indexStepLen, db.ChromList)
	for _, chrom := range db.ChromList {
		outfile := path.Join(transcriptDir, "chr"+chrom.Name+".idx.json")
		_, err := os.Stat(outfile)
		if os.IsExist(err) {
			continue
		}
		log.Printf("write %s", outfile)
		subTranscripts := transcripts.FilterByChrom(chrom.Name)
		subIndexes := indexes.FilterByChrom(chrom.Name)
		if subTranscripts.Len() > 0 && subIndexes.Len() > 0 {
			fo := utils.CreateFile(outfile)
			for _, _index := range subIndexes {
				_index.SetTranscript(subTranscripts)
				if len(_index.Transcripts) > 0 {
					contents := utils.ToJSON(_index)
					utils.WriteLine(fo, contents+"\n")
				}
			}
			utils.CloseFile(fo)
		}
	}
}
