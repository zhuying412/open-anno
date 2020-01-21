package core

import (
	"bufio"
	"bytes"
	"io"
	"os"
	"strings"
)

type Sequence []byte

func GetOne2Three(base byte) string {
	return AaOne2ThreeDict[base]
}

func (seq Sequence) Reverse() {
	length := len(seq)
	tmpStr := make(Sequence, length)
	copy(tmpStr, seq)
	for i, v := range tmpStr {
		seq[length-i-1] = v
	}
}

func (seq Sequence) GetLen() int {
	return len(seq)
}

func (seq Sequence) GetChar(index int) byte {
	return seq[index]
}

func (seq Sequence) GetSeq(index int, len int) Sequence {
	if len < 0 || index+len >= seq.GetLen() {
		return seq[index:]
	}
	return seq[index : index+len]
}

func (seq Sequence) GetIndex(base byte) int {
	return bytes.IndexByte(seq, base)
}

func (seq Sequence) IsStartswith(prefix Sequence) bool {
	return bytes.HasPrefix(seq, prefix)
}

func (seq Sequence) IsEndswith(suffix Sequence) bool {
	return bytes.HasSuffix(seq, suffix)
}

func (seq1 Sequence) IsEqual(seq2 Sequence) bool {
	return bytes.Equal(seq1, seq2)
}

func (seq Sequence) IsEmpty() bool {
	return len(seq) == 0
}

func (seq1 *Sequence) CopyFrom(seq2 Sequence) {
	*seq1 = make([]byte, seq2.GetLen())
	copy(*seq1, seq2)
}

func (seq *Sequence) RemoveOne(oldSeq Sequence) {
	*seq = bytes.Replace(*seq, oldSeq, Sequence(""), 1)
}

func (seq *Sequence) RemoveAll(oldSeq Sequence) {
	*seq = bytes.Replace(*seq, oldSeq, Sequence(""), -1)
}

func (seq *Sequence) Clear() {
	*seq = Sequence("")
}

func (seq *Sequence) PushChar(c byte) {
	*seq = append(*seq, c)
}

func (seq *Sequence) UnshiftChar(c byte) {
	tmp := make(Sequence, (*seq).GetLen()+1)
	tmp[0] = c
	copy(tmp[1:], *seq)
	*seq = tmp
}

func (seq1 *Sequence) Push(seq2 Sequence) {
	*seq1 = append(*seq1, seq2...)
}

func (seq1 *Sequence) Unshift(seq2 Sequence) {
	tmp := make(Sequence, (*seq1).GetLen()+seq2.GetLen())
	copy(tmp[0:seq2.GetLen()], seq2)
	copy(tmp[seq2.GetLen():], *seq1)
	*seq1 = tmp
}

func (seq Sequence) GetSnpSequence(pos int, alt byte) Sequence {
	newSeq := make(Sequence, seq.GetLen())
	copy(newSeq, seq)
	newSeq[pos-1] = alt
	return newSeq
}

func (seq Sequence) GetInsSequence(pos int, alt Sequence) Sequence {
	newSeq := make(Sequence, seq.GetLen()+alt.GetLen())
	copy(newSeq, seq.GetSeq(0, pos))
	copy(newSeq[pos:], alt)
	copy(newSeq[pos+len(alt):], seq.GetSeq(pos, -1))
	return newSeq
}

func (seq Sequence) GetDelSequence(lenL int, lenR int) Sequence {
	newSeq := make(Sequence, lenL+lenR)
	copy(newSeq, seq.GetSeq(0, lenL))
	copy(newSeq[lenL:], seq.GetSeq(seq.GetLen()-lenR, lenR))
	return newSeq
}

func (seq Sequence) Translate(isMt bool) Sequence {
	length := seq.GetLen()
	protein := make(Sequence, length/3)
	for i, j := 0, 0; i < length; i, j = i+3, j+1 {
		var codonDict map[string]byte
		if isMt {
			codonDict = CodonMtDict
		} else {
			codonDict = CodonDict
		}
		if aa, ok := codonDict[string(seq.GetSeq(i, 3))]; ok {
			protein[j] = aa
		}
	}
	return protein
}

func (protein Sequence) IsCmpl() bool {
	return bytes.IndexByte(protein, '*') != -1
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
