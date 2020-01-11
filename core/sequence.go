package core

import (
	"bytes"
)

type Sequence []byte

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
	if len < 0 {
		return seq[index:]
	}
	return seq[index : index+len]
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

func (seq1 *Sequence) Push(seq2 Sequence) {
	*seq1 = append(*seq1, seq2...)
}

func (seq1 *Sequence) Unshift(seq2 Sequence) {
	tmp := make(Sequence, (*seq1).GetLen()+seq2.GetLen())
	copy(tmp[0:seq2.GetLen()], seq2)
	copy(tmp[seq2.GetLen():], *seq1)
	*seq1 = tmp
}
