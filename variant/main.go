package variant

import (
	"OpenAnno/seq"
	"bufio"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

func NewVariant(chrom string, start string, end string, ref string, alt string) Variant {
	variant := Variant{Chrom: chrom, Ref: seq.Sequence(ref), Alt: seq.Sequence(alt)}
	var err error
	variant.Start, err = strconv.Atoi(start)
	if err != nil {
		log.Panic(err)
	}
	variant.End, err = strconv.Atoi(end)
	if err != nil {
		log.Panic(err)
	}
	return variant
}

func NewSnv(chrom string, pos int, ref seq.Sequence, alt seq.Sequence) Snv {
	if !ref.IsEmpty() || !alt.IsEmpty() && !ref.IsEqual(alt) {
		if ref.Startswith(alt) || ref.Endswith(alt) {
			if ref.Startswith(alt) {
				pos += alt.Len()
			}
			ref.Replace(alt, 1)
			alt.Clear()
		} else if alt.Startswith(ref) || alt.Endswith(ref) {
			if alt.Startswith(ref) {
				pos += ref.Len() - 1
			} else {
				pos += ref.Len() - alt.Len()
			}
			alt.Replace(ref, 1)
			ref.Clear()
		} else {
			var refRev, altRev seq.Sequence
			var subLen int
			refRev, altRev = ref, alt
			refRev.Reverse()
			altRev.Reverse()
			for i, subLen := 0, 0; i < ref.Len() && i < alt.Len(); i++ {
				if refRev.Base(i) != altRev.Base(i) {
					break
				}
				subLen++
			}
			ref = ref.SubSeq(0, ref.Len()-subLen)
			alt = alt.SubSeq(0, alt.Len()-subLen)
			for i, subLen := 0, 0; i < ref.Len() && i < alt.Len(); i++ {
				if ref.Base(i) != alt.Base(i) {
					break
				}
				subLen++
			}
			ref = ref.SubSeq(subLen, -1)
			alt = alt.SubSeq(subLen, -1)
			if subLen > 0 && ref.IsEmpty() {
				pos += subLen - 1
			} else {
				pos += subLen
			}
		}
	}
	snv := Snv{Variant: Variant{Chrom: chrom, Start: pos, End: pos, Ref: ref, Alt: alt}}
	snv.Chrom = chrom
	if snv.Chrom == "M" {
		snv.Chrom = "MT"
	}
	if snv.Ref.IsEmpty() {
		snv.End = snv.Start
		snv.Ref = "-"
	} else {
		snv.End = snv.Start + snv.Ref.Len() - 1
	}
	if snv.Alt.IsEmpty() {
		snv.Alt = "-"
	}
	return snv
}

func NewCnv(chrom string, start int, end int, copyNumber int) Cnv {
	var alt string
	if copyNumber > 1 {
		alt = "DUP"
	} else if copyNumber < 1 {
		alt = "DEL"
	} else {
		alt = "DIP"
	}
	return Cnv{
		Variant: Variant{
			Chrom: chrom,
			Start: start,
			End:   end,
			Ref:   "DIP",
			Alt:   seq.Sequence(alt),
		},
		CopyNumber: copyNumber,
	}
}

func ReadSnvFile(snvFile string) (Snvs, OtherInfoMap) {
	snvs := make(Snvs, 0)
	infoMap := make(OtherInfoMap)
	fi, err := os.Open(snvFile)
	if err == nil {
		log.Panic(err)
	}
	defer func(fp *os.File) {
		err := fp.Close()
		if err != nil {
			log.Panic(err)
		}
	}(fi)
	reader := bufio.NewReader(fi)
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Panic(err)
			}
		}
		line = strings.TrimSpace(line)
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		fields := strings.Split(line, "\t")
		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			log.Panic(err)
		}
		snv := NewSnv(fields[0], pos, seq.Sequence(fields[3]), seq.Sequence(fields[4]))
		snvs = append(snvs, snv)
		infoMap[snv.SN()] = NewOtherInfo(fields[5])
	}
	sort.Sort(snvs)
	return snvs, infoMap
}

func ReadCnvFile(cnvFile string) (Cnvs, OtherInfoMap) {
	cnvs := make(Cnvs, 0)
	infoMap := make(OtherInfoMap)
	fi, err := os.Open(cnvFile)
	if err == nil {
		log.Panic(err)
	}
	defer func(fp *os.File) {
		err := fp.Close()
		if err != nil {
			log.Panic(err)
		}
	}(fi)
	reader := bufio.NewReader(fi)
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Panic(err)
			}
		}
		line = strings.TrimSpace(line)
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		fields := strings.Split(line, "\t")
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			log.Panic(err)
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			log.Panic(err)
		}
		copyNumber, err := strconv.Atoi(fields[3])
		if err != nil {
			log.Panic(err)
		}
		cnv := NewCnv(fields[0], start, end, copyNumber)
		cnvs = append(cnvs, cnv)
		infoMap[cnv.SN()] = NewOtherInfo(fields[5])
	}
	sort.Sort(cnvs)
	return cnvs, infoMap
}
