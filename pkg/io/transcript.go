package io

import (
	"errors"
	"fmt"
	"log"
	"open-anno/pkg"
	"open-anno/pkg/schema"
	"open-anno/pkg/seq"
	"os/exec"
	"strconv"
	"strings"

	"github.com/brentp/faidx"
)

type TransScanner struct {
	Scanner[schema.Transcript]
}

func NewTransScanner(reader Reader) TransScanner {
	scanner := NewScanner[schema.Transcript](reader)
	return TransScanner{Scanner: scanner}
}

func (this TransScanner) Row() (schema.Transcript, error) {
	transcript := schema.Transcript{}
	fields := strings.Split(this.Text(), "\t")
	var name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, gene string
	var exonStarts, exonEnds []string
	switch len(fields) {
	case 16:
		name = fields[1]
		chrom = pkg.FormatChrom(fields[2])
		strand = fields[3]
		txStart = fields[4]
		txEnd = fields[5]
		cdsStart = fields[6]
		cdsEnd = fields[7]
		exonCount = fields[8]
		exonStarts = strings.Split(strings.Trim(fields[9], ","), ",")
		exonEnds = strings.Split(strings.Trim(fields[10], ","), ",")
		gene = fields[12]
	case 12:
		name = fields[0]
		chrom = pkg.FormatChrom(fields[1])
		strand = fields[2]
		txStart = fields[3]
		txEnd = fields[4]
		cdsStart = fields[5]
		cdsEnd = fields[6]
		exonCount = fields[7]
		exonStarts = strings.Split(strings.Trim(fields[8], ","), ",")
		exonEnds = strings.Split(strings.Trim(fields[9], ","), ",")
		gene = fields[0]
	default:
		return transcript, errors.New("unkown refGene format")
	}
	var err error
	transcript.Name = name
	transcript.Chrom = chrom
	transcript.Strand = strand
	transcript.Gene = gene
	transcript.TxStart, err = strconv.Atoi(txStart)
	if err != nil {
		return transcript, err
	}
	transcript.TxStart++
	transcript.TxEnd, err = strconv.Atoi(txEnd)
	if err != nil {
		return transcript, err
	}
	transcript.CdsStart, err = strconv.Atoi(cdsStart)
	if err != nil {
		return transcript, err
	}
	transcript.CdsStart++
	transcript.CdsEnd, err = strconv.Atoi(cdsEnd)
	if err != nil {
		return transcript, err
	}
	transcript.ExonCount, err = strconv.Atoi(exonCount)
	if err != nil {
		return transcript, err
	}
	transcript.ExonStarts = make([]int, transcript.ExonCount)
	transcript.ExonEnds = make([]int, transcript.ExonCount)
	for i := 0; i < transcript.ExonCount; i++ {
		transcript.ExonStarts[i], err = strconv.Atoi(exonStarts[i])
		if err != nil {
			return transcript, err
		}
		transcript.ExonStarts[i]++
		transcript.ExonEnds[i], err = strconv.Atoi(exonEnds[i])
		if err != nil {
			return transcript, err
		}
	}
	return transcript, err
}

func ReadGenePred(infile string) (schema.Transcripts, error) {
	transcripts := make(schema.Transcripts)
	fi, err := NewIoReader(infile)
	if err != nil {
		return transcripts, err
	}
	defer fi.Close()
	scanner := NewTransScanner(fi)
	for scanner.Scan() {
		transcript, err := scanner.Row()
		if err != nil {
			return transcripts, err
		}
		if _, ok := seq.GENOME[transcript.Chrom]; ok {
			transcripts[transcript.Name] = transcript
		}
	}
	return transcripts, err
}

func CreateAndIndexmRNA(transcripts schema.Transcripts, fai *faidx.Faidx, outfile string) error {
	writer, err := NewIoWriter(outfile)
	if err != err {
		return err
	}
	defer writer.Close()
	for _, trans := range transcripts {
		sequence, err := seq.Fetch(fai, trans.Chrom, trans.TxStart-1, trans.TxEnd)
		if err != nil {
			return err
		}
		sequence = strings.ToUpper(sequence)
		fmt.Fprintf(writer, ">%s:%s:%s\n%s\n", trans.Chrom, trans.Gene, trans.Name, sequence)
	}
	command := exec.Command("samtools", "faidx", outfile)
	err = command.Run()
	if err != nil {
		log.Print(err)
		log.Printf("Now you need run the command: 'samtools faidx %s'", outfile)
	}
	return nil
}
