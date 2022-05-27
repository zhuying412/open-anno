package gene

import (
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"open-anno/pkg/io/refgene"
	"open-anno/pkg/seq"
	"strings"
)

func getSubNcdna(cdna string, start int, end int, alt string, strand string, inSplcing bool) string {
	if inSplcing {
		nEnd := pkg.Min(end-start+1, len(alt))
		if strand == "+" {
			alt = alt[0:nEnd]
		} else {
			alt = seq.Reverse(seq.Reverse(alt)[0:nEnd])
		}
	}
	return seq.Substitute2(cdna, start, end, alt)
}

func setSubAAChange(anno SnvGeneBased, trans refgene.Transcript, cstart int, cend int, alt string, aashort bool) SnvGeneBased {
	cdna := trans.CDNA()
	ncdna := getSubNcdna(cdna, cstart, cend, alt, trans.Strand, strings.Contains(anno.Region, "splic"))
	if trans.Strand == "-" {
		alt = seq.RevComp(alt)
		cdna = seq.RevComp(cdna)
		ncdna = seq.RevComp(ncdna)
	}
	start := seq.DifferenceSimple(cdna, ncdna)
	if anno.NAChange == "" {
		if cstart == cend {
			anno.NAChange = fmt.Sprintf("c.%ddelins%s", start, alt)
		} else {
			anno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", start, start+cend-cstart, alt)
		}

	}
	protein := seq.Translate(cdna, trans.Chrom == "MT")
	nprotein := seq.Translate(ncdna, trans.Chrom == "MT")
	start, end1, end2 := seq.Difference(protein, nprotein)
	aa1 := protein[start-1 : end1]
	aa2 := nprotein[start-1 : end2]
	if (len(cdna)-len(ncdna))%3 == 0 {
		anno.Event = "sub_nonframeshift"
		if len(aa2) == 0 {
			if len(aa1) == 1 {
				anno.AAChange = fmt.Sprintf("p.%s%ddel", seq.AAName(aa1, aashort), start)
			} else {
				anno.AAChange = fmt.Sprintf(
					"p.%s%d_%s%ddel",
					seq.AAName(aa1[0], aashort),
					start,
					seq.AAName(aa1[len(aa1)-1], aashort),
					end1,
				)
			}
		} else {
			if len(aa1) == 1 {
				anno.AAChange = fmt.Sprintf("p.%s%ddelins%s", seq.AAName(aa1, aashort), start, seq.AAName(aa2, aashort))
			} else {
				anno.AAChange = fmt.Sprintf(
					"p.%s%d_%s%ddelins%s",
					seq.AAName(aa1[0], aashort),
					start,
					seq.AAName(aa1[len(aa1)-1], aashort),
					end1,
					seq.AAName(aa2, aashort))
			}
		}
	} else {
		if start < len(protein) {
			anno.Event = "sub_frameshift"
			var fs string
			fsi := strings.IndexByte(nprotein[start-1:], '*')
			if fsi == -1 {
				fs = "?"
			}
			if fsi != 0 {
				fs = fmt.Sprintf("%d", fsi+1)
			}
			anno.AAChange = fmt.Sprintf("p.%s%d%sfs*%s", seq.AAName(aa1[0], aashort), start, seq.AAName(aa2[0], aashort), fs)
		}
	}
	if protein[0] != nprotein[0] && protein[0] == 'M' {
		anno.Event += "_startloss"
	}
	if strings.Contains(anno.Region, "splic") {
		anno.Event += "_splicing"
	}
	return anno
}

func AnnoSub(snv io.Variant, trans refgene.Transcript, aashort bool) SnvGeneBased {
	cStart, cEnd, region1, region2, isExonSplicing := getDelCLen(trans, snv)
	cLen := trans.CLen()
	l := trans.CdsStart - pkg.Max(trans.TxStart, snv.Start)
	r := pkg.Min(trans.TxEnd, snv.End) - trans.CdsEnd
	anno := NewSnvGeneBased(trans, region1, region2)
	if isExonSplicing {
		anno.Region = "exonic_splicing"
	}
	if snv.Start < trans.CdsStart && snv.End > trans.CdsEnd {
		// snv包含了整个编码区，即整个编码区被删除
		// ...+++,,,+++...
		//   |---------|
		//|---------------|
		if trans.Strand == "+" {
			// anno.NAChange = fmt.Sprintf("c.-%d_+%ddel%s", trans.CdsStart-snv.Start, snv.End-trans.CdsEnd, snv.Ref)
			if l == 0 {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", 1, cLen, snv.Alt)
				} else {
					anno.NAChange = fmt.Sprintf("c.%d_+%ddelins%s", 1, r, snv.Alt)
				}
			} else {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.-%d_%ddelins%s", l, cLen, snv.Alt)
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_+%ddelins%s", l, r, snv.Alt)
				}
			}
		} else {
			// anno.NAChange = fmt.Sprintf("c.1_+%ddel%s", snv.End-trans.CdsEnd, trans.CdsStart-snv.Start, seq.RevComp(snv.Ref))
			if l == 0 {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", 1, cLen, seq.RevComp(snv.Alt))
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_%ddelins%s", r, cLen, seq.RevComp(snv.Alt))
				}
			} else {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_+%ddelins%s", 1, l, seq.RevComp(snv.Alt))
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_+%ddelins%s", r, l, seq.RevComp(snv.Alt))
				}
			}
		}
		anno.Region = "transcript"
		anno.Region2 = "transcript"
		anno.Event = "CNV"
	} else if snv.Start < trans.CdsStart && snv.End < trans.CdsStart {
		// snv发生在整个编码区左边
		// |-|...+++,,,+++...
		//   |-|
		if trans.Strand == "+" {
			anno.NAChange = fmt.Sprintf("c.-%d_-%ddelins%s", l, trans.CdsStart-snv.End, snv.Alt)
		} else {
			anno.NAChange = fmt.Sprintf("c.+%d_+%ddelins%s", trans.CdsStart-snv.End, l, seq.RevComp(snv.Alt))
		}
	} else if snv.Start > trans.CdsEnd && snv.End > trans.CdsEnd {
		// snv发生在整个编码区左边
		// ...+++,,,+++...|-|
		//              |-|
		if trans.Strand == "+" {
			anno.NAChange = fmt.Sprintf("c.+%d_+%ddelins%s", snv.Start-trans.CdsEnd, r, snv.Alt)
		} else {
			anno.NAChange = fmt.Sprintf("c.-%d_-%ddelins%s", snv.End-trans.CdsEnd, r, seq.RevComp(snv.Alt))
		}
	} else if trans.CdsStart > snv.Start && trans.CdsStart <= snv.End && trans.CdsEnd >= snv.End {
		if region2.Type == refgene.RType_CDS {
			// ...+++,,,+++...
			//  |--|
			//|----|
			if trans.Strand == "+" {
				if l == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", 1, cEnd, snv.Alt)
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_%ddelins%s", l, cEnd, snv.Alt)
				}
			} else {
				if l == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", cLen-cEnd+1, cLen, seq.RevComp(snv.Alt))
				} else {
					anno.NAChange = fmt.Sprintf("c.%d_+%ddelins%s", cLen-cEnd+1, l, seq.RevComp(snv.Alt))
				}
			}
		} else {
			// intron
			// ...,,,+++...
			//  |--|
			//|----|
			if trans.Strand == "+" {
				if l == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%d+%ddelins%s", 1, cEnd, snv.End-region2.Start+1, snv.Alt)
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_%d+%ddelins%s", l, cEnd, snv.End-region2.Start+1, snv.Alt)
				}
			} else {
				if l == 0 {
					anno.NAChange = fmt.Sprintf("c.%d-%d_%ddelins%s", cLen-cEnd+1, snv.End-region2.Start+1, cLen, seq.RevComp(snv.Alt))
				} else {
					anno.NAChange = fmt.Sprintf("c.%d-%d_+%ddelins%s", cLen-cEnd+1, snv.End-region2.Start+1, l, seq.RevComp(snv.Alt))
				}
			}
		}
		anno = setSubAAChange(anno, trans, 1, cEnd, snv.Alt, aashort)
	} else if trans.CdsStart > snv.Start && trans.CdsStart <= snv.End && trans.CdsEnd >= snv.End {
		if region1.Type == refgene.RType_CDS {
			// ...,,,+++...
			//        |--|
			//        |----|
			if trans.Strand == "+" {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", cStart, cLen, snv.Alt)
				} else {
					anno.NAChange = fmt.Sprintf("c.%d_+%ddelins%s", cStart, r, snv.Alt)
				}
			} else {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", 1, cLen-cStart+1, seq.RevComp(snv.Alt))
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_%ddelins%s", r, cLen-cStart+1, seq.RevComp(snv.Alt))
				}

			}
			anno = setSubAAChange(anno, trans, cStart, cLen, snv.Alt, aashort)
		} else {
			// intron
			// ...+++,,,...
			//        |--|
			//        |----|
			if trans.Strand == "+" {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d-%d_%ddelins%s", cStart+1, region1.End-snv.Start+1, cLen, snv.Alt)
				} else {
					anno.NAChange = fmt.Sprintf("c.%d-%d_+%ddelins%s", cStart+1, region1.End-snv.Start+1, r, snv.Alt)
				}
			} else {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%d+%ddelins%s", 1, cLen-cStart, region1.End-snv.Start+1, seq.RevComp(snv.Alt))
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_%d+%ddelins%s", r, cLen-cStart, region1.End-snv.Start+1, seq.RevComp(snv.Alt))
				}

			}
			anno = setSubAAChange(anno, trans, cStart+1, cLen, snv.Alt, aashort)
		}

	} else {
		if region1.Equal(region2) {
			if region1.Type == refgene.RType_CDS {
				// ...++++++,,,...
				//     |--|
				anno = setSubAAChange(anno, trans, cStart, cEnd, snv.Alt, aashort)
			} else {
				// ...+++,,,,,,...
				//        |--|
				dist1s, dist1e := snv.Start-region1.Start+1, snv.End-region1.Start+1
				dist2s, dist2e := region1.End-region1.Start+1, region1.End-snv.End+1
				dist1, dist2 := dist1s, dist2e
				if pkg.Min(dist1, dist2) <= 2 {
					anno.Event = "splicing"
				}
				if trans.Strand == "+" {
					if dist1 <= dist2 {
						anno.NAChange = fmt.Sprintf("c.%d+%d_%d+%ddelins%s", cStart, dist1s, cStart, dist1e, snv.Alt)
					} else {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%d-%ddelins%s", cStart+1, dist2s, cStart+1, dist2e, snv.Alt)
					}
				} else {
					if dist1 <= dist2 {
						nclen := cLen - cStart + 1
						anno.NAChange = fmt.Sprintf("c.%d-%d_%d-%ddelins%s", nclen, dist1e, nclen, dist1s, seq.RevComp(snv.Alt))
					} else {
						nclen := cLen - cStart
						anno.NAChange = fmt.Sprintf("c.%d+%d_%d+%ddelins%s", nclen, dist1e, nclen, dist1s, seq.RevComp(snv.Alt))
					}
				}
			}
		} else {
			if region1.Type == refgene.RType_CDS {
				if region2.Type != refgene.RType_CDS {
					// ...+++,,,+++...
					//     |--|
					if trans.Strand == "+" {
						anno.NAChange = fmt.Sprintf("c.%d_%d+%ddelins%s", cStart, cEnd, snv.End-region2.Start+1, snv.Alt)
					} else {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%ddelins%s", cLen-cEnd+1, snv.End-region2.Start+1, cLen-cStart+1, seq.RevComp(snv.Alt))
					}
				}
				anno = setSubAAChange(anno, trans, cStart, cEnd, snv.Alt, aashort)
			} else {
				if region2.Type == refgene.RType_CDS {
					// ...+++,,,+++...
					//        |--|
					if trans.Strand == "+" {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%ddelins%s", cStart+1, region1.End-snv.Start+1, cEnd, snv.Alt)
					} else {
						anno.NAChange = fmt.Sprintf("c.%d_%d+%ddelins%s", cLen-cEnd+1, cLen-cStart, region1.End-snv.Start+1, seq.RevComp(snv.Alt))
					}
				} else {
					// ...+++,,,+++,,,+++...
					//        |-----|
					if trans.Strand == "+" {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%d+%ddelins%s", cStart+1, region1.End-snv.Start+1, cEnd, snv.End-region2.Start+1, snv.Alt)
					} else {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%d+%ddelins%s", cLen-cEnd+1, snv.End-region2.Start+1, cLen-cStart, region1.End-snv.Start+1, seq.RevComp(snv.Alt))
					}
				}
				anno = setSubAAChange(anno, trans, cStart+1, cEnd, snv.Alt, aashort)
			}
		}
	}
	return anno
}
