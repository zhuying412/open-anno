package seq

import (
	"bytes"
	"log"
	"reflect"
	"strings"
)

type Base = byte

type Sequence string

func (s Sequence) Len() int {
	return len(s)
}

func (s Sequence) Base(i int) Base {
	return s[i]
}

func (s Sequence) SubSeq(index int, len int) (sub Sequence) {
	if len < 0 || index+len >= s.Len() {
		len = s.Len() - index
	}
	return s[index : index+len]
}

func (s Sequence) Find(sub interface{}) (i int) {
	switch sub.(type) {
	case string:
		i = strings.Index(string(s), sub.(string))
	case int32:
		i = strings.Index(string(s), string([]byte{byte(sub.(int32))}))
	case Base:
		i = strings.Index(string(s), string(sub.(Base)))
	case Sequence:
		i = strings.Index(string(s), string(sub.(Sequence)))
	default:
		log.Panicf("%v(%v) is not string, Base or Sequence", reflect.TypeOf(sub).Name(), sub)
	}
	return i
}

func (s Sequence) Startswith(sub Sequence) bool {
	return strings.HasPrefix(string(s), string(sub))
}

func (s Sequence) Endswith(sub Sequence) bool {
	return strings.HasSuffix(string(s), string(sub))
}

func (s Sequence) IsEqual(t Sequence) bool {
	return string(s) == string(t)
}

func (s Sequence) IsEmpty() bool {
	return s.Len() == 0
}

func (s Sequence) Translate(isMt bool) Sequence {
	var buffer bytes.Buffer
	for i := 0; i < s.Len(); i += 3 {
		codon := codonMap
		if isMt {
			codon = codonMtMap
		}
		if aa, ok := codon[s.SubSeq(i, 3)]; ok {
			buffer.WriteByte(aa)
		}
	}
	return Sequence(buffer.String())
}

func (s Sequence) ChangeWithSnp(pos int, alt Base) Sequence {
	var buffer bytes.Buffer
	buffer.WriteString(string(s.SubSeq(0, pos-1)))
	buffer.WriteByte(alt)
	buffer.WriteString(string(s.SubSeq(pos, -1)))
	return Sequence(buffer.String())
}

func (s Sequence) ChangeWithIns(pos int, alt Sequence) Sequence {
	var buffer bytes.Buffer
	buffer.WriteString(string(s.SubSeq(0, pos)))
	buffer.WriteString(string(alt))
	buffer.WriteString(string(s.SubSeq(pos, -1)))
	return Sequence(buffer.String())
}

func (s Sequence) ChangeWithDel(lenL int, lenR int) Sequence {
	var buffer bytes.Buffer
	buffer.WriteString(string(s.SubSeq(0, lenL)))
	buffer.WriteString(string(s.SubSeq(s.Len()-lenR, lenR)))
	return Sequence(buffer.String())
}

func (s Sequence) ProteinOne2Tree() string {
	var buffer bytes.Buffer
	for i := 0; i < s.Len(); i++ {
		buffer.WriteString(AAMap[s.Base(i)])
	}
	return buffer.String()
}

func (s *Sequence) Reverse() {
	var buffer bytes.Buffer
	for i := s.Len() - 1; i >= 0; i-- {
		buffer.WriteByte(s.Base(i))
	}
	*s = Sequence(buffer.String())
}

func (s *Sequence) Replace(sub interface{}, count int) {
	switch sub.(type) {
	case string:
		*s = Sequence(strings.Replace(string(*s), sub.(string), "", count))
	case int32:
		*s = Sequence(strings.Replace(string(*s), string([]byte{byte(sub.(int32))}), "", count))
	case Base:
		*s = Sequence(strings.Replace(string(*s), string(sub.(Base)), "", count))
	case Sequence:
		*s = Sequence(strings.Replace(string(*s), string(sub.(Sequence)), "", count))
	default:
		log.Panicf("%v(%v) is not string, Base or Sequence", reflect.TypeOf(sub).Name(), sub)
	}
}

func (s *Sequence) Clear() {
	*s = ""
}

func (s *Sequence) Join(seqs ...Sequence) {
	var buffer bytes.Buffer
	buffer.WriteString(string(*s))
	for _, seq := range seqs {
		buffer.WriteString(string(seq))
	}
	*s = Sequence(buffer.String())
}

func (s Sequence) IsProteinCmpl() bool {
	return s.Find('*') != -1
}
