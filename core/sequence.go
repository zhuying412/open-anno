package core

import (
	"bufio"
	"bytes"
	"io"
	"os"
	"strings"
)

type Base = byte
type Sequence string

func GetOne2Three(base Base) string {
	return AaOne2ThreeDict[base]
}

func (seq Sequence) String() string {
	return string(seq)
}

func (seq *Sequence) Reverse() {
	var buffer bytes.Buffer
	for i := seq.GetLen() - 1; i >= 0; i-- {
		buffer.WriteByte(seq.GetChar(i))
	}
	*seq = Sequence(buffer.String())
}

func (seq Sequence) GetLen() int {
	return len(seq)
}

func (seq Sequence) GetChar(index int) Base {
	return seq[index]
}

func (seq Sequence) GetSeq(index int, len int) Sequence {
	if len < 0 || index+len >= seq.GetLen() {
		return seq[index:]
	}
	return seq[index : index+len]
}

func (seq Sequence) GetIndex(base Base) int {
	return strings.IndexByte(seq.String(), base)
}

func (seq Sequence) IsStartswith(prefix Sequence) bool {
	return strings.HasPrefix(seq.String(), prefix.String())
}

func (seq Sequence) IsEndswith(suffix Sequence) bool {
	return strings.HasSuffix(seq.String(), suffix.String())
}

func (seq1 Sequence) IsEqual(seq2 Sequence) bool {
	return seq1 == seq2
}

func (seq Sequence) IsEmpty() bool {
	return seq == ""
}

func (seq *Sequence) RemoveOne(substr Sequence) {
	*seq = Sequence(strings.Replace(seq.String(), substr.String(), "", 1))
}

func (seq *Sequence) RemoveAll(substr Sequence) {
	*seq = Sequence(strings.Replace(seq.String(), substr.String(), "", -1))
}

func (seq *Sequence) Clear() {
	*seq = ""
}

func (seq *Sequence) Join(sequences []Sequence) {
	var buffer bytes.Buffer
	for _, sequence := range sequences {
		buffer.WriteString(sequence.String())
	}
	*seq = Sequence(buffer.String())
}

func (seq Sequence) GetSnpSequence(pos int, alt Base) Sequence {
	var buffer bytes.Buffer
	buffer.WriteString(seq.GetSeq(0, pos-1).String())
	buffer.WriteByte(alt)
	buffer.WriteString(seq.GetSeq(pos, -1).String())
	return Sequence(buffer.String())
}

func (seq Sequence) GetInsSequence(pos int, alt Sequence) Sequence {
	var buffer bytes.Buffer
	buffer.WriteString(seq.GetSeq(0, pos).String())
	buffer.WriteString(alt.String())
	buffer.WriteString(seq.GetSeq(pos, -1).String())
	return Sequence(buffer.String())
}

func (seq Sequence) GetDelSequence(lenL int, lenR int) Sequence {
	var buffer bytes.Buffer
	buffer.WriteString(seq.GetSeq(0, lenL).String())
	buffer.WriteString(seq.GetSeq(seq.GetLen()-lenR, lenR).String())
	return Sequence(buffer.String())
}

func (seq Sequence) Translate(isMt bool) Sequence {
	var buffer bytes.Buffer
	for i := 0; i < seq.GetLen(); i += 3 {
		var codonDict map[string]byte
		if isMt {
			codonDict = CodonMtDict
		} else {
			codonDict = CodonDict
		}
		if aa, ok := codonDict[seq.GetSeq(i, 3).String()]; ok {
			buffer.WriteByte(aa)
		}
	}
	return Sequence(buffer.String())
}

func (protein Sequence) IsCmpl() bool {
	return strings.IndexByte(protein.String(), '*') != -1
}

func (protein Sequence) GetOne2Tree() string {
	var buffer bytes.Buffer
	for i := 0; i < protein.GetLen(); i++ {
		buffer.WriteString(GetOne2Three(protein[i]))
	}
	return buffer.String()
}

type Fasta map[string]Sequence

func (fasta Fasta) Read(fastaFile string) {
	if fp, err := os.Open(fastaFile); err == nil {
		defer fp.Close()
		reader := bufio.NewReader(fp)
		var name, seq bytes.Buffer
		for {
			if line, err := reader.ReadBytes('\n'); err == nil {
				line = bytes.TrimSpace(line)
				if len(line) == 0 {
					continue
				}
				if line[0] == '>' {
					if name.Len() != 0 {
						_name := strings.Split(name.String(), " ")[0]
						fasta[_name] = Sequence(seq.String())
					}
					name.Reset()
					seq.Reset()
					name.Write(line[1:])
				} else {
					seq.Write(line)
				}

			} else {
				if err == io.EOF {
					break
				} else {
					panic(err.Error())
				}
			}
		}
		if name.Len() != 0 {
			_name := strings.Split(name.String(), " ")[0]
			fasta[_name] = Sequence(seq.String())
		}
	} else {
		panic(err.Error())
	}
}
