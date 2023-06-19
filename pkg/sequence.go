package pkg

import (
	"bytes"
	"fmt"
	"regexp"
)

// RevComp 反向互补
func Reverse(sequence string) string {
	var buffer bytes.Buffer
	for i := len(sequence) - 1; i >= 0; i-- {
		buffer.WriteByte(sequence[i])
	}
	return buffer.String()
}

// RevComp 反向互补
func RevComp(sequence string) string {
	var buffer bytes.Buffer
	for i := len(sequence) - 1; i >= 0; i-- {
		buffer.WriteByte(ATGCs[sequence[i]])
	}
	return buffer.String()
}

// Translate 翻译蛋白
func Translate(sequence string, mt bool) string {
	var buffer bytes.Buffer
	length := len(sequence)
	for i := 0; i+3 <= length; i += 3 {
		na := sequence[i : i+3]
		table := NAtoAATable
		if mt {
			table = MitoNAtoAATable
		}
		if aa, ok := table[na]; ok {
			buffer.WriteByte(aa)
		} else {
			buffer.WriteByte('X')
		}
	}
	if length%3 > 0 {
		buffer.WriteByte('X')
	}
	return buffer.String()
}

// AAName 氨基酸名称
func AAName[T byte | string](bases T, aashort bool) string {
	var buffer bytes.Buffer
	sequence := string(bases)
	for i := 0; i < len(sequence); i++ {
		if aashort {
			buffer.WriteByte(sequence[i])
		} else {
			buffer.WriteString(AAMap[sequence[i]])
		}
	}
	return buffer.String()
}

// Substitute 替换碱基
func Substitute(sequence string, pos int, base string) string {
	var buffer bytes.Buffer
	buffer.WriteString(sequence[0 : pos-1])
	buffer.WriteString(base)
	buffer.WriteString(string(sequence[pos:]))
	return buffer.String()
}

// Insert 插入碱基
func Insert(sequence string, pos int, bases string) string {
	var buffer bytes.Buffer
	buffer.WriteString(sequence[0:pos])
	buffer.WriteString(bases)
	buffer.WriteString(sequence[pos:])
	return buffer.String()
}

// Delete 删除碱基
func Delete(sequence string, start int, end int) string {
	var buffer bytes.Buffer
	buffer.WriteString(sequence[0 : start-1])
	buffer.WriteString(sequence[end:])
	return buffer.String()
}

// Substitute2 替换碱基
func Substitute2(sequence string, start int, end int, alt string) string {
	var buffer bytes.Buffer
	buffer.WriteString(sequence[0 : start-1])
	buffer.WriteString(alt)
	buffer.WriteString(sequence[end:])
	return buffer.String()
}

// DupUnit Dup的单位元件：如 ATGATGATG->ATG
func DupUnit(sequence string) string {
	var unit string
	for i := 0; i < len(sequence); i++ {
		unit = sequence[0 : i+1]
		match, _ := regexp.MatchString(fmt.Sprintf(`^(%s)+$`, unit), sequence)
		if match {
			break
		}
	}
	return unit
}

// DifferenceSimple 比较两个序列，返回第一个差异位置
func DifferenceSimple(sequence1 string, sequence2 string) int {
	minLen := Min(len(sequence1), len(sequence2))
	pos := 0
	for i := 0; i < minLen && sequence1[i] == sequence2[i]; i++ {
		pos++
	}
	return pos + 1
}

// Difference 比较两个序列，返回差异序列的第一个差异位置、及各自最后一个差异位置
func Difference(sequence1 string, sequence2 string) (int, int, int) {
	minLen := Min(len(sequence1), len(sequence2))
	lLen, rLen := 0, 0
	for i := 0; i < minLen && sequence1[i] == sequence2[i]; i++ {
		lLen++
	}
	for i, j := len(sequence1)-1, len(sequence2)-1; i >= 0 && j >= 0 && sequence1[i] == sequence2[j]; i, j = i-1, j-1 {
		rLen++
	}
	if lLen+rLen > minLen {
		rLen = minLen - lLen
	}
	return lLen + 1, len(sequence1) - rLen, len(sequence2) - rLen
	// return lLen, rLen
}
