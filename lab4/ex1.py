"""
Translate a coding-region DNA/RNA sequence into an amino acid sequence
"""
from typing import Tuple
import sys



CODON_TABLE = {
	# Phenylalanine
	'UUU': 'F', 'UUC': 'F',
	# Leucine
	'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
	# Isoleucine
	'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
	# Methionine (Start)
	'AUG': 'M',
	# Valine
	'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
	# Serine
	'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',
	# Proline
	'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
	# Threonine
	'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
	# Alanine
	'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	# Tyrosine
	'UAU': 'Y', 'UAC': 'Y',
	# Histidine
	'CAU': 'H', 'CAC': 'H',
	# Glutamine
	'CAA': 'Q', 'CAG': 'Q',
	# Asparagine
	'AAU': 'N', 'AAC': 'N',
	# Lysine
	'AAA': 'K', 'AAG': 'K',
	# Aspartic acid
	'GAU': 'D', 'GAC': 'D',
	# Glutamic acid
	'GAA': 'E', 'GAG': 'E',
	# Cysteine
	'UGU': 'C', 'UGC': 'C',
	# Tryptophan
	'UGG': 'W',
	# Arginine
	'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
	# Glycine
	'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

STOP_CODONS = {'UAA', 'UAG', 'UGA'}


def _normalize_seq(seq: str) -> str:
	
	s = seq.strip().upper()
	
	s = s.replace('T', 'U')
	
	s = ''.join([c for c in s if c in 'AUCG'])
	return s


def translate_coding_region(seq: str, find_first_start: bool = True) -> Tuple[str, int]:
	
	rna = _normalize_seq(seq)
	if find_first_start:
		start = rna.find('AUG')
		if start == -1:
			
			start = 0
	else:
		start = 0

	protein_chars = []
	i = start
	while i + 3 <= len(rna):
		codon = rna[i:i+3]
		if codon in STOP_CODONS:
			break
		aa = CODON_TABLE.get(codon, 'X') 
		protein_chars.append(aa)
		i += 3

	return ''.join(protein_chars), start


def _main(argv):
	if len(argv) >= 2:
		seq = argv[1]
	else:
		
		seq = sys.stdin.read().strip()

	if not seq:
		print('Usage: python ex1.py <sequence>  (or provide sequence via stdin)')
		return 1

	protein, start = translate_coding_region(seq, find_first_start=True)
	print(protein)
	return 0


if __name__ == '__main__':
	raise SystemExit(_main(sys.argv))

