package db

import (
	"OpenAnno/pkg/utils"
	"fmt"
	"log"
	"os"
	"strings"
)

func SplitDatabase(databaseFile string, outdir string) {
	utils.CreateDir(outdir)
	log.Printf("read %s", databaseFile)
	handlerMap := make(map[string]*os.File)
	fi, reader := utils.OpenFile(databaseFile)
	defer utils.CloseFile(fi)
	var header string
	for {
		line, isEof := utils.ReadLine(reader, 0)
		if isEof {
			break
		}
		if line[0] == '#' {
			if header == "" {
				header = strings.Trim(line, "#")
			}
			continue
		}
		fields := strings.Split(line, "\t")
		chrom := fields[0]
		if _, ok := handlerMap[chrom]; !ok {
			handlerMap[chrom] = utils.CreateFile(fmt.Sprintf("%s/chr%s.txt", outdir, chrom))
			utils.WriteLine(handlerMap[chrom], "#"+header+"\n")
		}
		utils.WriteLine(handlerMap[chrom], line+"\n")
	}
	for _, handler := range handlerMap {
		utils.CloseFile(handler)
	}
}
