#! /bin/bash

infile=$1
outfile=$(echo $infile|sed 's:.txt$:.sorted.txt:g')

sort -k 1,1 -k 2n,2 -k 3n,3 -k 4,4 -k 5,5 -o $outfile $infile && rm $infile