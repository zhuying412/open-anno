package core

var (
	ChromList = []string{
		"1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
		"11", "12", "13", "14", "15", "16", "17", "18",
		"19", "20", "21", "22", "X", "Y", "MT",
	}
	ChromOrderDict = map[string]int{
		"1": 1, "2": 2, "3": 3, "4": 4, "5": 5, "6": 6, "7": 7, "8": 8, "9": 9, "10": 10,
		"11": 11, "12": 12, "13": 13, "14": 14, "15": 15, "16": 16, "17": 17, "18": 18,
		"19": 19, "20": 20, "21": 21, "22": 22, "X": 23, "Y": 24, "MT": 25,
	}
	ChromLengthDict = map[string]int{
		"1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260, "6": 171115067,
		"7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
		"13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753, "17": 81195210, "18": 78077248,
		"19": 59128983, "20": 63025520, "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
		"MT": 16569,
	}
	CodonDicT = map[string]byte{
		"TTT": 'F', "TTC": 'F', "TTA": 'L', "TTG": 'L', "TCT": 'S', "TCC": 'S', "TCA": 'S', "TCG": 'S',
		"TAT": 'Y', "TAC": 'Y', "TAA": '*', "TAG": '*', "TGT": 'C', "TGC": 'C', "TGA": '*', "TGG": 'W',
		"CTT": 'L', "CTC": 'L', "CTA": 'L', "CTG": 'L', "CCT": 'P', "CCC": 'P', "CCA": 'P', "CCG": 'P',
		"CAT": 'H', "CAC": 'H', "CAA": 'Q', "CAG": 'Q', "CGT": 'R', "CGC": 'R', "CGA": 'R', "CGG": 'R',
		"ATT": 'I', "ATC": 'I', "ATA": 'I', "ATG": 'M', "ACT": 'T', "ACC": 'T', "ACA": 'T', "ACG": 'T',
		"AAT": 'N', "AAC": 'N', "AAA": 'K', "AAG": 'K', "AGT": 'S', "AGC": 'S', "AGA": 'R', "AGG": 'R',
		"GTT": 'V', "GTC": 'V', "GTA": 'V', "GTG": 'V', "GCT": 'A', "GCC": 'A', "GCA": 'A', "GCG": 'A',
		"GAT": 'D', "GAC": 'D', "GAA": 'E', "GAG": 'E', "GGT": 'G', "GGC": 'G', "GGA": 'G', "GGG": 'G',
	}
	CodonMtDicT = map[string]byte{
		"TTT": 'F', "TTC": 'F', "TTA": 'L', "TTG": 'L', "TCT": 'S', "TCC": 'S', "TCA": 'S', "TCG": 'S',
		"TAT": 'Y', "TAC": 'Y', "TAA": '*', "TAG": '*', "TGT": 'C', "TGC": 'C', "TGA": 'W', "TGG": 'W',
		"CTT": 'L', "CTC": 'L', "CTA": 'L', "CTG": 'L', "CCT": 'P', "CCC": 'P', "CCA": 'P', "CCG": 'P',
		"CAT": 'H', "CAC": 'H', "CAA": 'Q', "CAG": 'Q', "CGT": 'R', "CGC": 'R', "CGA": 'R', "CGG": 'R',
		"ATT": 'I', "ATC": 'I', "ATA": 'M', "ATG": 'M', "ACT": 'T', "ACC": 'T', "ACA": 'T', "ACG": 'T',
		"AAT": 'N', "AAC": 'N', "AAA": 'K', "AAG": 'K', "AGT": 'S', "AGC": 'S', "AGA": '*', "AGG": '*',
		"GTT": 'V', "GTC": 'V', "GTA": 'V', "GTG": 'V', "GCT": 'A', "GCC": 'A', "GCA": 'A', "GCG": 'A',
		"GAT": 'D', "GAC": 'D', "GAA": 'E', "GAG": 'E', "GGT": 'G', "GGC": 'G', "GGA": 'G', "GGG": 'G',
	}
	AaOne2ThreeDict = map[byte]string{
		'A': "Ala", 'R': "Arg", 'N': "Asn", 'D': "Asp", 'C': "Cys", 'Q': "Gln", 'E': "Glu",
		'G': "Gly", 'H': "His", 'I': "Ile", 'L': "Leu", 'K': "Lys", 'M': "Met", 'F': "Phe",
		'P': "Pro", 'S': "Ser", 'T': "Thr", 'W': "Trp", 'Y': "Tyr", 'V': "Val", 'X': "Ter",
		'*': "Ter",
	}
)
