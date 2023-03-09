package gene

import (
	"fmt"
	"open-anno/anno"
	"open-anno/pkg"
	"strings"
)

func getDelCLen(trans pkg.Transcript, snv anno.AnnoVariant) (int, int, pkg.Region, pkg.Region, bool) {
	var cStart, cEnd int
	var region1, region2 pkg.Region
	var hasCDS, hasIntron bool
	for _, region := range trans.Regions {
		if region.End < snv.Start {
			if region.Type == pkg.RType_CDS {
				cStart += region.End - region.Start + 1
			}
		}
		if region.Start <= snv.Start && snv.Start <= region.End {
			if region.Type == pkg.RType_CDS {
				cStart += snv.Start - region.Start + 1
			}
			region1 = region
		}
		if region.End < snv.End {
			if region.Type == pkg.RType_CDS {
				cEnd += region.End - region.Start + 1
			}
		}
		if region.Start <= snv.End && snv.End <= region.End {
			if region.Type == pkg.RType_CDS {
				cEnd += snv.End - region.Start + 1
			}
			region2 = region
		}
		if region.Start <= snv.End && region.End >= snv.Start {
			if region.Type == pkg.RType_CDS {
				hasCDS = true
			}
			if region.Type == pkg.RType_INTRON {
				hasIntron = true
			}
		}
	}

	return cStart, cEnd, region1, region2, hasCDS && hasIntron
}

func setDelAAChange(transAnno TransAnno, trans pkg.Transcript, cstart int, cend int) TransAnno {
	cdna := trans.CDNA()
	ncdna := pkg.Delete(cdna, cstart, cend)
	if trans.Strand == "-" {
		cdna = pkg.RevComp(cdna)
		ncdna = pkg.RevComp(ncdna)
	}
	start := pkg.DifferenceSimple(cdna, ncdna)
	alt := cdna[start-1 : start+cend-cstart]
	if transAnno.NAChange == "" {
		if cstart == cend {
			transAnno.NAChange = fmt.Sprintf("c.%ddel%s", start, alt)
		} else {
			transAnno.NAChange = fmt.Sprintf("c.%d_%ddel%s", start, start+cend-cstart, alt)
		}

	}
	protein := pkg.Translate(cdna, trans.Chrom == "MT")
	nprotein := pkg.Translate(ncdna, trans.Chrom == "MT")
	start, end1, end2 := pkg.Difference(protein, nprotein)
	aa1 := protein[start-1 : end1]
	aa2 := nprotein[start-1 : end2]
	if (len(cdna)-len(ncdna))%3 == 0 {
		transAnno.Event = "del_inframe"
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
			if aa1[len(aa1)-1] == '*' {
				transAnno.AAChange += "ext*?"
				transAnno.Event += "_stoploss"
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
			transAnno.Event = "del_frameshift"
			if len(aa2) == 0 {
				transAnno.AAChange = fmt.Sprintf("p.%s%dfs", pkg.AAName(aa1[0], AA_SHORT), start)
			} else {
				if aa2[0] == '*' {
					transAnno.AAChange = fmt.Sprintf("p.%s%d*", pkg.AAName(aa1[0], AA_SHORT), start)
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
	}
	if protein[0] == 'M' && nprotein[0] != 'M' {
		transAnno.Event += "_startloss"
	}
	if strings.Contains(transAnno.Region, "splic") {
		transAnno.Event += "_splicing"
	}
	if protein[len(protein)-1] == '*' && strings.IndexByte(nprotein, '*') == -1 {
		transAnno.Event += "_stoploss"
	}
	return transAnno
}

func AnnoDel(snv anno.AnnoVariant, trans pkg.Transcript) TransAnno {
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
					transAnno.NAChange = fmt.Sprintf("c.%d_%ddel", 1, cLen)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d_+%ddel", 1, r)
				}
			} else {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.-%d_%ddel", l, cLen)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_+%ddel", l, r)
				}
			}
		} else {
			// transAnno.NAChange = fmt.Sprintf("c.1_+%ddel%s", snv.End-trans.CdsEnd, trans.CdsStart-snv.Start, pkg.RevComp(snv.Ref))
			if l == 0 {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%ddel", 1, cLen)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_%ddel", r, cLen)
				}
			} else {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_+%ddel", 1, l)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_+%ddel", r, l)
				}
			}
		}
		transAnno.AAChange = fmt.Sprintf("c.%d_%ddel", 1, cLen)
		transAnno.Region = "transcript"
		transAnno.Region2 = "transcript"
		transAnno.Event = "CNV"
	} else if snv.Start < trans.CdsStart && snv.End < trans.CdsStart {
		// snv发生在整个编码区左边
		// |-|...+++,,,+++...
		//   |-|
		ll := pkg.Abs(trans.CdsStart-snv.End-l) + 1
		if trans.Strand == "+" {
			transAnno.NAChange = fmt.Sprintf("c.-%d_-%ddel%s", l, trans.CdsStart-snv.End, snv.Ref[len(snv.Ref)-ll:])
		} else {
			transAnno.NAChange = fmt.Sprintf("c.+%d_+%ddel%s", trans.CdsStart-snv.End, l, pkg.RevComp(snv.Ref[len(snv.Ref)-ll:]))
		}
	} else if snv.Start > trans.CdsEnd && snv.End > trans.CdsEnd {
		// snv发生在整个编码区左边
		// ...+++,,,+++...|-|
		//              |-|
		ll := pkg.Abs(snv.Start-trans.CdsEnd-r) + 1
		if trans.Strand == "+" {
			transAnno.NAChange = fmt.Sprintf("c.+%d_+%ddel%s", snv.Start-trans.CdsEnd, r, snv.Ref[0:ll])
		} else {
			transAnno.NAChange = fmt.Sprintf("c.-%d_-%ddel%s", snv.End-trans.CdsEnd, r, pkg.RevComp(snv.Ref[0:ll]))
		}
	} else if trans.CdsStart > snv.Start && trans.CdsStart <= snv.End && trans.CdsEnd >= snv.End {
		ll := snv.End - trans.CdsStart + 1 + l
		if region2.Type == pkg.RType_CDS {
			// ...+++,,,+++...
			//  |--|
			//|----|
			if trans.Strand == "+" {
				if l == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%ddel%s", 1, cEnd, snv.Ref[len(snv.Ref)-ll:])
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_%ddel%s", l, cEnd, snv.Ref[len(snv.Ref)-ll:])
				}
			} else {

				if l == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%ddel%s", cLen-cEnd+1, cLen, pkg.RevComp(snv.Ref[len(snv.Ref)-ll:]))
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d_+%ddel%s", cLen-cEnd+1, l, pkg.RevComp(snv.Ref[len(snv.Ref)-ll:]))
				}
			}
		} else {
			// intron
			// ...,,,+++...
			//  |--|
			//|----|
			if trans.Strand == "+" {
				if l == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%d+%ddel%s", 1, cEnd, snv.End-region2.Start+1, snv.Ref)
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_%d+%ddel%s", l, cEnd, snv.End-region2.Start+1, snv.Ref)
				}
			} else {
				if l == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d-%d_%ddel%s", cLen-cEnd+1, snv.End-region2.Start+1, cLen, pkg.RevComp(snv.Ref))
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d-%d_+%ddel%s", cLen-cEnd+1, snv.End-region2.Start+1, l, pkg.RevComp(snv.Ref))
				}
			}
		}
		transAnno = setDelAAChange(transAnno, trans, 1, cEnd)
	} else if trans.CdsStart > snv.Start && trans.CdsStart <= snv.End && trans.CdsEnd >= snv.End {
		ll := trans.CdsEnd - snv.End + 1 + r
		if region1.Type == pkg.RType_CDS {
			// ...,,,+++...
			//        |--|
			//        |----|
			if trans.Strand == "+" {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%ddel%s", cStart, cLen, snv.Ref[0:ll])
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d_+%ddel%s", cStart, r, snv.Ref[0:ll])
				}
			} else {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%ddel%s", 1, cLen-cStart+1, pkg.RevComp(snv.Ref[0:ll]))
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_%ddel%s", r, cLen-cStart+1, pkg.RevComp(snv.Ref[0:ll]))
				}

			}
			transAnno = setDelAAChange(transAnno, trans, cStart, cLen)
		} else {
			// intron
			// ...+++,,,...
			//        |--|
			//        |----|
			if trans.Strand == "+" {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d-%d_%ddel%s", cStart+1, region1.End-snv.Start+1, cLen, snv.Ref[0:ll])
				} else {
					transAnno.NAChange = fmt.Sprintf("c.%d-%d_+%ddel%s", cStart+1, region1.End-snv.Start+1, r, snv.Ref[0:ll])
				}
			} else {
				if r == 0 {
					transAnno.NAChange = fmt.Sprintf("c.%d_%d+%ddel%s", 1, cLen-cStart, region1.End-snv.Start+1, pkg.RevComp(snv.Ref[0:ll]))
				} else {
					transAnno.NAChange = fmt.Sprintf("c.-%d_%d+%ddel%s", r, cLen-cStart, region1.End-snv.Start+1, pkg.RevComp(snv.Ref[0:ll]))
				}

			}
			transAnno = setDelAAChange(transAnno, trans, cStart+1, cLen)
		}

	} else {
		if region1.Equal(region2) {
			if region1.Type == pkg.RType_CDS {
				// ...++++++,,,...
				//     |--|
				transAnno = setDelAAChange(transAnno, trans, cStart, cEnd)
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
						transAnno.NAChange = fmt.Sprintf("c.%d+%d_%d+%ddel%s", cStart, dist1s, cStart, dist1e, snv.Ref)
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%d-%ddel%s", cStart+1, dist2s, cStart+1, dist2e, snv.Ref)
					}
				} else {
					if dist1 <= dist2 {
						nclen := cLen - cStart + 1
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%d-%ddel%s", nclen, dist1e, nclen, dist1s, pkg.RevComp(snv.Ref))
					} else {
						nclen := cLen - cStart
						transAnno.NAChange = fmt.Sprintf("c.%d+%d_%d+%ddel%s", nclen, dist1e, nclen, dist1s, pkg.RevComp(snv.Ref))
					}
				}
			}
		} else {
			if region1.Type == pkg.RType_CDS {
				if region2.Type != pkg.RType_CDS {
					// ...+++,,,+++...
					//     |--|
					if trans.Strand == "+" {
						transAnno.NAChange = fmt.Sprintf("c.%d_%d+%ddel%s", cStart, cEnd, snv.End-region2.Start+1, snv.Ref)
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%ddel%s", cLen-cEnd+1, snv.End-region2.Start+1, cLen-cStart+1, pkg.RevComp(snv.Ref))
					}
				}
				transAnno = setDelAAChange(transAnno, trans, cStart, cEnd)
			} else {
				if region2.Type == pkg.RType_CDS {
					// ...+++,,,+++...
					//        |--|
					if trans.Strand == "+" {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%ddel%s", cStart+1, region1.End-snv.Start+1, cEnd, snv.Ref)
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%d_%d+%ddel%s", cLen-cEnd+1, cLen-cStart, region1.End-snv.Start+1, pkg.RevComp(snv.Ref))
					}
				} else {
					// ...+++,,,+++,,,+++...
					//        |-----|
					if trans.Strand == "+" {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%d+%ddel%s", cStart+1, region1.End-snv.Start+1, cEnd, snv.End-region2.Start+1, snv.Ref)
					} else {
						transAnno.NAChange = fmt.Sprintf("c.%d-%d_%d+%ddel%s", cLen-cEnd+1, snv.End-region2.Start+1, cLen-cStart, region1.End-snv.Start+1, pkg.RevComp(snv.Ref))
					}
				}
				transAnno = setDelAAChange(transAnno, trans, cStart+1, cEnd)
			}
		}
	}
	return transAnno
}
