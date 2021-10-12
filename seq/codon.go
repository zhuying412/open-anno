package seq

var codonMap = map[Sequence]Base{
	"TTT": 'F', "TTC": 'F', "TTA": 'L', "TTG": 'L', "TCT": 'S', "TCC": 'S', "TCA": 'S', "TCG": 'S',
	"TAT": 'Y', "TAC": 'Y', "TAA": '*', "TAG": '*', "TGT": 'C', "TGC": 'C', "TGA": '*', "TGG": 'W',
	"CTT": 'L', "CTC": 'L', "CTA": 'L', "CTG": 'L', "CCT": 'P', "CCC": 'P', "CCA": 'P', "CCG": 'P',
	"CAT": 'H', "CAC": 'H', "CAA": 'Q', "CAG": 'Q', "CGT": 'R', "CGC": 'R', "CGA": 'R', "CGG": 'R',
	"ATT": 'I', "ATC": 'I', "ATA": 'I', "ATG": 'M', "ACT": 'T', "ACC": 'T', "ACA": 'T', "ACG": 'T',
	"AAT": 'N', "AAC": 'N', "AAA": 'K', "AAG": 'K', "AGT": 'S', "AGC": 'S', "AGA": 'R', "AGG": 'R',
	"GTT": 'V', "GTC": 'V', "GTA": 'V', "GTG": 'V', "GCT": 'A', "GCC": 'A', "GCA": 'A', "GCG": 'A',
	"GAT": 'D', "GAC": 'D', "GAA": 'E', "GAG": 'E', "GGT": 'G', "GGC": 'G', "GGA": 'G', "GGG": 'G',
}

var codonMtMap = map[Sequence]Base{
	"TTT": 'F', "TTC": 'F', "TTA": 'L', "TTG": 'L', "TCT": 'S', "TCC": 'S', "TCA": 'S', "TCG": 'S',
	"TAT": 'Y', "TAC": 'Y', "TAA": '*', "TAG": '*', "TGT": 'C', "TGC": 'C', "TGA": 'W', "TGG": 'W',
	"CTT": 'L', "CTC": 'L', "CTA": 'L', "CTG": 'L', "CCT": 'P', "CCC": 'P', "CCA": 'P', "CCG": 'P',
	"CAT": 'H', "CAC": 'H', "CAA": 'Q', "CAG": 'Q', "CGT": 'R', "CGC": 'R', "CGA": 'R', "CGG": 'R',
	"ATT": 'I', "ATC": 'I', "ATA": 'M', "ATG": 'M', "ACT": 'T', "ACC": 'T', "ACA": 'T', "ACG": 'T',
	"AAT": 'N', "AAC": 'N', "AAA": 'K', "AAG": 'K', "AGT": 'S', "AGC": 'S', "AGA": '*', "AGG": '*',
	"GTT": 'V', "GTC": 'V', "GTA": 'V', "GTG": 'V', "GCT": 'A', "GCC": 'A', "GCA": 'A', "GCG": 'A',
	"GAT": 'D', "GAC": 'D', "GAA": 'E', "GAG": 'E', "GGT": 'G', "GGC": 'G', "GGA": 'G', "GGG": 'G',
}

var AAMap = map[Base]string{
	'A': "Ala", 'R': "Arg", 'N': "Asn", 'D': "Asp", 'C': "Cys", 'Q': "Gln", 'E': "Glu",
	'G': "Gly", 'H': "His", 'I': "Ile", 'L': "Leu", 'K': "Lys", 'M': "Met", 'F': "Phe",
	'P': "Pro", 'S': "Ser", 'T': "Thr", 'W': "Trp", 'Y': "Tyr", 'V': "Val", 'X': "Ter",
	'*': "Ter",
}
