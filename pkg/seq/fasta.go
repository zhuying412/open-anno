package seq

import (
	"OpenAnno/pkg/utils"
	"bytes"
	"fmt"
	"strings"
)

type Fasta map[string]Sequence

func (f *Fasta) Set(name string, seq Sequence) {
	name = strings.Split(name, " ")[0]
	(*f)[name] = seq
}

func ReadFastaFile(path string) (fasta Fasta) {
	fasta = make(Fasta)
	fi, reader := utils.OpenFile(path)
	defer utils.CloseFile(fi)
	var name, sequence bytes.Buffer
	for {
		line, isEof := utils.ReadLine(reader, 0)
		if isEof {
			break
		}
		if line[0] == '>' {
			if name.Len() != 0 {
				fasta.Set(name.String(), Sequence(sequence.String()))
			}
			name.Reset()
			sequence.Reset()
			name.WriteString(line[1:])
		} else {
			sequence.WriteString(strings.ToUpper(line))
		}
	}
	if name.Len() != 0 {
		fasta.Set(name.String(), Sequence(sequence.String()))
	}
	return fasta
}

func CreateFastaFile(fasta Fasta, path string) {
	fo := utils.CreateFile(path)
	defer utils.CloseFile(fo)
	for sn, sequence := range fasta {
		utils.WriteLine(fo, fmt.Sprintf(">%s\n%s\n", sn, sequence))
	}
}
