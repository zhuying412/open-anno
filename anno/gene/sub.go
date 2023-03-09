package gene

import (
	"fmt"
	"open-anno/anno"
	"open-anno/pkg"
	"strings"
)

func getSubNcdna(cdna string, start int, end int, alt string, strand string, inSplcing bool) string {
	if inSplcing {
		nEnd := pkg.Min(end-start+1, len(alt))
		if strand == "+" {
			alt = alt[0:nEnd]
		} else {
			alt = pkg.Reverse(pkg.Reverse(alt)[0:nEnd])
		}
	}
	return pkg.Substitute2(cdna, start, end, alt)
}

func setSubAAChange(transAnno TransAnno, trans pkg.Transcript, cstart int, cend int, alt string) TransAnno {
	cdna := trans.CDNA()
	ncdna := getSubNcdna(cdna, cstart, cend, alt, trans.Strand, strings.Contains(transAnno.Region, "splic"))
	if trans.Strand == "-" {
		alt = pkg.RevComp(alt)
		cdna = pkg.RevComp(cdna)
		ncdna = pkg.RevComp(ncdna)
	}
	start := pkg.DifferenceSimple(cdna, ncdna)
	if transAnno.NAChange == "" {
		if cstart == cend {
			transAnno.NAChange = fmt.Sprintf("c.%ddelins%s", start, alt)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", start, start+cend-cstart, alt)
		}

	}
	protein := pkg.Translate(cdna, trans.Chrom == "MT")
	nprotein := pkg.Translate(ncdna, trans.Chrom == "MT")
	start, end1, end2 := pkg.Difference(protein, nprotein)
	aa1 := protein[start-1 : end1]
	aa2 := nprotein[start-1 : end2]
	if (len(cdna)-len(ncdna))%3 == 0 {
		transAnno.Event = "sub_inframe"
		if len(aa2) == 0 {
			if len(aa1) == 1 {
				transAnno.AAChange = fmt.Sprintf("p.%s%ddel", pkg.AAName(aa1, AA_SHORT), start)
			} else {
				transAnno.AAChange = fmt.Sprintf(
					"p.%s%d_%s%ddel",
					pkg.AAName(aa1[0], AA_SHORT),
					start,
					pkg.AAName(aa1[len(aa1)-1], AA_SHORT),
					end1,
				)
			}
		} else {
			if len(aa1) == 1 {
				transAnno.AAChange = fmt.Sprintf("p.%s%ddelins%s", pkg.AAName(aa1, AA_SHORT), start, pkg.AAName(aa2, AA_SHORT))
			} else {
				transAnno.AAChange = fmt.Sprintf(
					"p.%s%d_%s%ddelins%s",
					pkg.AAName(aa1[0], AA_SHORT),
					start,
					pkg.AAName(aa1[len(aa1)-1], AA_SHORT),
					end1,
					pkg.AAName(aa2, AA_SHORT))
			}
		}
	} else {
		if start < len(protein) {
			transAnno.Event = "sub_frameshift"
			if aa2[0] == '*' {
				transAnno.AAChange = fmt.Sprintf("p.%s%dfs", pkg.AAName(aa1[0], AA_SHORT), start)
			} else {
				var fs string
				fsi := strings.IndexByte(nprotein[start-1:], '*')
				if fsi == -1 {
					fs = "?"
				}
				if fsi != 0 {
					fs = fmt.Sprintf("%d", fsi+1)
				}

				transAnno.AAChange = fmt.Sprintf("p.%s%d%sfs*%s", pkg.AAName(aa1[0], AA_SHORT), start, pkg.AAName(aa2[0], AA_SHORT), fs)
			}
		}
	}
	if protein[0] != nprotein[0] && protein[0] == 'M' {
		transAnno.Event += "_startloss"
	}
	if strings.Contains(transAnno.Region, "splic") {
		transAnno.Event += "_splicing"
	}
	return transAnno
}

func AnnoSub(snv anno.AnnoVariant, trans pkg.Transcript) TransAnno {
	cStart, cEnd, region1, region2, isExonSplicing := getDelCLen(trans, snv)
	cLen := trans.CLen()
	l := trans.CdsStart - pkg.Max(trans.TxStart, snv.Start)
	r := pkg.Min(trans.TxEnd, snv.End) - trans.CdsEnd
	transAnno := NewTransAnno(trans, region1, region2)
	if isExonSplicing {
		transAnno.Region = "exonic_splicing"
	}
	if snv.Start < trans.CdsStart && snv.End > trans.CdsEnd {
		// snv包含了整个编码区，即整个编码区被删除
		// ...+++,,,+++...
		//   |---------|
		//|---------------|
		if trans.Strand == "+" {
			// transAnno.NAChange = fmt.Sprintf("c.-%d_+%ddel%s", trans.CdsStart-snv.Start, snv.End-trans.CdsEnd, snv.Ref)
			if l == 0 {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", 1, cLen, snv.Alt)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d_+%ddelins%s", 1, r, snv.Alt)
				}
			} else {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.-%d_%ddelins%s", l, cLen, snv.Alt)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_+%ddelins%s", l, r, snv.Alt)
				}
			}
		} else {
			// transAnno.NAChange = fmt.Sprintf("c.1_+%ddel%s", snv.End-trans.CdsEnd, trans.CdsStart-snv.Start, pkg.RevComp(snv.Ref))
			if l == 0 {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", 1, cLen, pkg.RevComp(snv.Alt))
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_%ddelins%s", r, cLen, pkg.RevComp(snv.Alt))
				}
			} else {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_+%ddelins%s", 1, l, pkg.RevComp(snv.Alt))
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_+%ddelins%s", r, l, pkg.RevComp(snv.Alt))
				}
			}
		}
		transAnno.Region = "transcript"
		transAnno.Region2 = "transcript"
		transAnno.Event = "CNV"
	} else if snv.Start < trans.CdsStart && snv.End < trans.CdsStart {
		// snv发生在整个编码区左边
		// |-|...+++,,,+++...
		//   |-|
		if trans.Strand == "+" {
			transAnno.NAChange = fmt.Sprintf("c.-%d_-%ddelins%s", l, trans.CdsStart-snv.End, snv.Alt)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.+%d_+%ddelins%s", trans.CdsStart-snv.End, l, pkg.RevComp(snv.Alt))
		}
	} else if snv.Start > trans.CdsEnd && snv.End > trans.CdsEnd {
		// snv发生在整个编码区左边
		// ...+++,,,+++...|-|
		//              |-|
		if trans.Strand == "+" {
			transAnno.NAChange = fmt.Sprintf("c.+%d_+%ddelins%s", snv.Start-trans.CdsEnd, r, snv.Alt)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.-%d_-%ddelins%s", snv.End-trans.CdsEnd, r, pkg.RevComp(snv.Alt))
		}
	} else if trans.CdsStart > snv.Start && trans.CdsStart <= snv.End && trans.CdsEnd >= snv.End {
		if region2.Type == pkg.RType_CDS {
			// ...+++,,,+++...
			//  |--|
			//|----|
			if trans.Strand == "+" {
				if l == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", 1, cEnd, snv.Alt)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_%ddelins%s", l, cEnd, snv.Alt)
				}
			} else {
				if l == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", cLen-cEnd+1, cLen, pkg.RevComp(snv.Alt))
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d_+%ddelins%s", cLen-cEnd+1, l, pkg.RevComp(snv.Alt))
				}
			}
		} else {
			// intron
			// ...,,,+++...
			//  |--|
			//|----|
			if trans.Strand == "+" {
				if l == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%d+%ddelins%s", 1, cEnd, snv.End-region2.Start+1, snv.Alt)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_%d+%ddelins%s", l, cEnd, snv.End-region2.Start+1, snv.Alt)
				}
			} else {
				if l == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d-%d_%ddelins%s", cLen-cEnd+1, snv.End-region2.Start+1, cLen, pkg.RevComp(snv.Alt))
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d-%d_+%ddelins%s", cLen-cEnd+1, snv.End-region2.Start+1, l, pkg.RevComp(snv.Alt))
				}
			}
		}
		transAnno = setSubAAChange(transAnno, trans, 1, cEnd, snv.Alt)
	} else if trans.CdsStart > snv.Start && trans.CdsStart <= snv.End && trans.CdsEnd >= snv.End {
		if region1.Type == pkg.RType_CDS {
			// ...,,,+++...
			//        |--|
			//        |----|
			if trans.Strand == "+" {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", cStart, cLen, snv.Alt)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d_+%ddelins%s", cStart, r, snv.Alt)
				}
			} else {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%ddelins%s", 1, cLen-cStart+1, pkg.RevComp(snv.Alt))
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_%ddelins%s", r, cLen-cStart+1, pkg.RevComp(snv.Alt))
				}

			}
			transAnno = setSubAAChange(transAnno, trans, cStart, cLen, snv.Alt)
		} else {
			// intron
			// ...+++,,,...
			//        |--|
			//        |----|
			if trans.Strand == "+" {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d-%d_%ddelins%s", cStart+1, region1.End-snv.Start+1, cLen, snv.Alt)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d-%d_+%ddelins%s", cStart+1, region1.End-snv.Start+1, r, snv.Alt)
				}
			} else {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%d+%ddelins%s", 1, cLen-cStart, region1.End-snv.Start+1, pkg.RevComp(snv.Alt))
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_%d+%ddelins%s", r, cLen-cStart, region1.End-snv.Start+1, pkg.RevComp(snv.Alt))
				}

			}
			transAnno = setSubAAChange(transAnno, trans, cStart+1, cLen, snv.Alt)
		}

	} else {
		if region1.Equal(region2) {
			if region1.Type == pkg.RType_CDS {
				// ...++++++,,,...
				//     |--|
				transAnno = setSubAAChange(transAnno, trans, cStart, cEnd, snv.Alt)
			} else {
				// ...+++,,,,,,...
				//        |--|
				dist1s, dist1e := snv.Start-region1.Start+1, snv.End-region1.Start+1
				dist2s, dist2e := region1.End-region1.Start+1, region1.End-snv.End+1
				dist1, dist2 := dist1s, dist2e
				if pkg.Min(dist1, dist2) <= 2 {
					transAnno.Event = "splicing"
					transAnno.Region = "splicing"
				}
				if trans.Strand == "+" {
					if dist1 <= dist2 {
						transAnno.NAChange = fmt.Sprintf("c.%d+%d_%d+%ddelins%s", cStart, dist1s, cStart, dist1e, snv.Alt)
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%d-%ddelins%s", cStart+1, dist2s, cStart+1, dist2e, snv.Alt)
					}
				} else {
					if dist1 <= dist2 {
						nclen := cLen - cStart + 1
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%d-%ddelins%s", nclen, dist1e, nclen, dist1s, pkg.RevComp(snv.Alt))
					} else {
						nclen := cLen - cStart
						transAnno.NAChange = fmt.Sprintf("c.%d+%d_%d+%ddelins%s", nclen, dist1e, nclen, dist1s, pkg.RevComp(snv.Alt))
					}
				}
			}
		} else {
			if region1.Type == pkg.RType_CDS {
				if region2.Type != pkg.RType_CDS {
					// ...+++,,,+++...
					//     |--|
					if trans.Strand == "+" {
						transAnno.NAChange = fmt.Sprintf("c.%d_%d+%ddelins%s", cStart, cEnd, snv.End-region2.Start+1, snv.Alt)
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%ddelins%s", cLen-cEnd+1, snv.End-region2.Start+1, cLen-cStart+1, pkg.RevComp(snv.Alt))
					}
				}
				transAnno = setSubAAChange(transAnno, trans, cStart, cEnd, snv.Alt)
			} else {
				if region2.Type == pkg.RType_CDS {
					// ...+++,,,+++...
					//        |--|
					if trans.Strand == "+" {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%ddelins%s", cStart+1, region1.End-snv.Start+1, cEnd, snv.Alt)
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%d_%d+%ddelins%s", cLen-cEnd+1, cLen-cStart, region1.End-snv.Start+1, pkg.RevComp(snv.Alt))
					}
				} else {
					// ...+++,,,+++,,,+++...
					//        |-----|
					if trans.Strand == "+" {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%d+%ddelins%s", cStart+1, region1.End-snv.Start+1, cEnd, snv.End-region2.Start+1, snv.Alt)
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%d+%ddelins%s", cLen-cEnd+1, snv.End-region2.Start+1, cLen-cStart, region1.End-snv.Start+1, pkg.RevComp(snv.Alt))
					}
				}
				transAnno = setSubAAChange(transAnno, trans, cStart+1, cEnd, snv.Alt)
			}
		}
	}
	return transAnno
}
