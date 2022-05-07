package gene

import (
	"fmt"
	"open-anno/pkg"
	"open-anno/pkg/io"
	"open-anno/pkg/io/refgene"
	"open-anno/pkg/seq"
	"strings"
)

func getCLen(trans refgene.Transcript, snv io.Variant) (int, int, refgene.Region, refgene.Region, bool) {
	var cStart, cEnd int
	var region1, region2 refgene.Region
	var hasCDS, hasIntron bool
	for _, region := range trans.Regions {
		if region.End < snv.Start {
			if region.Type == refgene.RType_CDS {
				cStart += region.End - region.Start + 1
			}
		}
		if region.Start <= snv.Start && snv.Start <= region.End {
			if region.Type == refgene.RType_CDS {
				cStart += snv.Start - region.Start + 1
			}
			region1 = region
		}
		if region.End < snv.End {
			if region.Type == refgene.RType_CDS {
				cEnd += region.End - region.Start + 1
			}
		}
		if region.Start <= snv.End && snv.End <= region.End {
			if region.Type == refgene.RType_CDS {
				cEnd += snv.End - region.Start + 1
			}
			region2 = region
		}
		if region.Start <= snv.End && region.End >= snv.Start {
			if region.Type == refgene.RType_CDS {
				hasCDS = true
			}
			if region.Type == refgene.RType_INTRON {
				hasIntron = true
			}
		}
	}

	return cStart, cEnd, region1, region2, hasCDS && hasIntron
}

func setAAChange(anno SnvGeneBased, trans refgene.Transcript, cstart int, cend int, aashort bool) SnvGeneBased {
	cdna := trans.CDNA()
	ncdna := seq.Delete(cdna, cstart, cend)
	start := seq.DifferenceSimple(cdna, ncdna)
	alt := cdna[start-1 : start+cend-cstart]
	if anno.NAChange == "" {
		if cstart == cend {
			anno.NAChange = fmt.Sprintf("c.%ddel%s", start, alt)
		} else {
			anno.NAChange = fmt.Sprintf("c.%d_%ddel%s", start, start+cend-cstart, alt)
		}

	}
	if trans.Strand == "-" {
		cdna = seq.RevComp(cdna)
		ncdna = seq.RevComp(ncdna)
	}
	protein := seq.Translate(cdna, trans.Chrom == "MT")
	nprotein := seq.Translate(ncdna, trans.Chrom == "MT")
	start, end1, end2 := seq.Difference(protein, nprotein)
	aa1 := protein[start-1 : end1]
	aa2 := nprotein[start-1 : end2]
	if (len(cdna)-len(ncdna))%3 == 0 {
		anno.Event = "del_nonframeshift"
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
			anno.Event = "del_frameshift"
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

func AnnoDel(snv io.Variant, trans refgene.Transcript, aashort bool) SnvGeneBased {
	cStart, cEnd, region1, region2, isExonSplicing := getCLen(trans, snv)
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
					anno.NAChange = fmt.Sprintf("c.%d_%ddel", 1, cLen)
				} else {
					anno.NAChange = fmt.Sprintf("c.%d_+%ddel", 1, r)
				}
			} else {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.-%d_%ddel", l, cLen)
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_+%ddel", l, r)
				}
			}
		} else {
			// anno.NAChange = fmt.Sprintf("c.1_+%ddel%s", snv.End-trans.CdsEnd, trans.CdsStart-snv.Start, seq.RevComp(snv.Ref))
			if l == 0 {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%ddel", 1, cLen)
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_%ddel", r, cLen)
				}
			} else {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_+%ddel", 1, l)
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_+%ddel", r, l)
				}
			}
		}
		anno.AAChange = fmt.Sprintf("c.%d_%ddel", 1, cLen)
		anno.Region = "transcript"
		anno.Region2 = "transcript"
		anno.Event = "CNV"
	} else if snv.Start < trans.CdsStart && snv.End < trans.CdsStart {
		// snv发生在整个编码区左边
		// |-|...+++,,,+++...
		//   |-|
		ll := pkg.Abs(trans.CdsStart-snv.End-l) + 1
		if trans.Strand == "+" {
			anno.NAChange = fmt.Sprintf("c.-%d_-%ddel%s", l, trans.CdsStart-snv.End, snv.Ref[len(snv.Ref)-ll:])
		} else {
			anno.NAChange = fmt.Sprintf("c.+%d_+%ddel%s", trans.CdsStart-snv.End, l, seq.RevComp(snv.Ref[len(snv.Ref)-ll:]))
		}
	} else if snv.Start > trans.CdsEnd && snv.End > trans.CdsEnd {
		// snv发生在整个编码区左边
		// ...+++,,,+++...|-|
		//              |-|
		ll := pkg.Abs(snv.Start-trans.CdsEnd-r) + 1
		if trans.Strand == "+" {
			anno.NAChange = fmt.Sprintf("c.+%d_+%ddel%s", snv.Start-trans.CdsEnd, r, snv.Ref[0:ll])
		} else {
			anno.NAChange = fmt.Sprintf("c.-%d_-%ddel%s", snv.End-trans.CdsEnd, r, seq.RevComp(snv.Ref[0:ll]))
		}
	} else if trans.CdsStart > snv.Start && trans.CdsStart <= snv.End && trans.CdsEnd >= snv.End {
		ll := snv.End - trans.CdsStart + 1 + l
		if region2.Type == refgene.RType_CDS {
			// ...+++,,,+++...
			//  |--|
			//|----|
			if trans.Strand == "+" {
				if l == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%ddel%s", 1, cEnd, snv.Ref[len(snv.Ref)-ll:])
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_%ddel%s", l, cEnd, snv.Ref[len(snv.Ref)-ll:])
				}
			} else {

				if l == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%ddel%s", cLen-cEnd+1, cLen, seq.RevComp(snv.Ref[len(snv.Ref)-ll:]))
				} else {
					anno.NAChange = fmt.Sprintf("c.%d_+%ddel%s", cLen-cEnd+1, l, seq.RevComp(snv.Ref[len(snv.Ref)-ll:]))
				}
			}
		} else {
			// intron
			// ...,,,+++...
			//  |--|
			//|----|
			if trans.Strand == "+" {
				if l == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%d+%ddel%s", 1, cEnd, snv.End-region2.Start+1, snv.Ref)
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_%d+%ddel%s", l, cEnd, snv.End-region2.Start+1, snv.Ref)
				}
			} else {
				if l == 0 {
					anno.NAChange = fmt.Sprintf("c.%d-%d_%ddel%s", cLen-cEnd+1, snv.End-region2.Start+1, cLen, seq.RevComp(snv.Ref))
				} else {
					anno.NAChange = fmt.Sprintf("c.%d-%d_+%ddel%s", cLen-cEnd+1, snv.End-region2.Start+1, l, seq.RevComp(snv.Ref))
				}
			}
		}
		anno = setAAChange(anno, trans, 1, cEnd, aashort)
	} else if trans.CdsStart > snv.Start && trans.CdsStart <= snv.End && trans.CdsEnd >= snv.End {
		ll := trans.CdsEnd - snv.End + 1 + r
		if region1.Type == refgene.RType_CDS {
			// ...,,,+++...
			//        |--|
			//        |----|
			if trans.Strand == "+" {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%ddel%s", cStart, cLen, snv.Ref[0:ll])
				} else {
					anno.NAChange = fmt.Sprintf("c.%d_+%ddel%s", cStart, r, snv.Ref[0:ll])
				}
			} else {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%ddel%s", 1, cLen-cStart+1, seq.RevComp(snv.Ref[0:ll]))
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_%ddel%s", r, cLen-cStart+1, seq.RevComp(snv.Ref[0:ll]))
				}

			}
			anno = setAAChange(anno, trans, cStart, cLen, aashort)
		} else {
			// intron
			// ...+++,,,...
			//        |--|
			//        |----|
			if trans.Strand == "+" {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d-%d_%ddel%s", cStart+1, region1.End-snv.Start+1, cLen, snv.Ref[0:ll])
				} else {
					anno.NAChange = fmt.Sprintf("c.%d-%d_+%ddel%s", cStart+1, region1.End-snv.Start+1, r, snv.Ref[0:ll])
				}
			} else {
				if r == 0 {
					anno.NAChange = fmt.Sprintf("c.%d_%d+%ddel%s", 1, cLen-cStart, region1.End-snv.Start+1, seq.RevComp(snv.Ref[0:ll]))
				} else {
					anno.NAChange = fmt.Sprintf("c.-%d_%d+%ddel%s", r, cLen-cStart, region1.End-snv.Start+1, seq.RevComp(snv.Ref[0:ll]))
				}

			}
			anno = setAAChange(anno, trans, cStart+1, cLen, aashort)
		}

	} else {
		if region1.Equal(region2) {
			if region1.Type == refgene.RType_CDS {
				// ...++++++,,,...
				//     |--|
				anno = setAAChange(anno, trans, cStart, cEnd, aashort)
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
						anno.NAChange = fmt.Sprintf("c.%d+%d_%d+%ddel%s", cStart, dist1s, cStart, dist1e, snv.Ref)
					} else {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%d-%ddel%s", cStart+1, dist2s, cStart+1, dist2e, snv.Ref)
					}
				} else {
					if dist1 <= dist2 {
						nclen := cLen - cStart + 1
						anno.NAChange = fmt.Sprintf("c.%d-%d_%d-%ddel%s", nclen, dist1e, nclen, dist1s, seq.RevComp(snv.Ref))
					} else {
						nclen := cLen - cStart
						anno.NAChange = fmt.Sprintf("c.%d+%d_%d+%ddel%s", nclen, dist1e, nclen, dist1s, seq.RevComp(snv.Ref))
					}
				}
			}
		} else {
			if region1.Type == refgene.RType_CDS {
				if region2.Type != refgene.RType_CDS {
					// ...+++,,,+++...
					//     |--|
					if trans.Strand == "+" {
						anno.NAChange = fmt.Sprintf("c.%d_%d+%ddel%s", cStart, cEnd, snv.End-region2.Start+1, snv.Ref)
					} else {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%ddel%s", cLen-cEnd+1, snv.End-region2.Start+1, cLen-cStart+1, seq.RevComp(snv.Ref))
					}
				}
				anno = setAAChange(anno, trans, cStart, cEnd, aashort)
			} else {
				if region2.Type == refgene.RType_CDS {
					// ...+++,,,+++...
					//        |--|
					if trans.Strand == "+" {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%ddel%s", cStart+1, region1.End-snv.Start+1, cEnd, snv.Ref)
					} else {
						anno.NAChange = fmt.Sprintf("c.%d_%d+%ddel%s", cLen-cEnd+1, cLen-cStart, region1.End-snv.Start+1, seq.RevComp(snv.Ref))
					}
				} else {
					// ...+++,,,+++,,,+++...
					//        |-----|
					if trans.Strand == "+" {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%d+%ddel%s", cStart+1, region1.End-snv.Start+1, cEnd, snv.End-region2.Start+1, snv.Ref)
					} else {
						anno.NAChange = fmt.Sprintf("c.%d-%d_%d+%ddel%s", cLen-cEnd+1, snv.End-region2.Start+1, cLen-cStart, region1.End-snv.Start+1, seq.RevComp(snv.Ref))
					}
				}
				anno = setAAChange(anno, trans, cStart+1, cEnd, aashort)
			}
		}
	}
	return anno
}
