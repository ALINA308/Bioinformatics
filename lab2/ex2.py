def find_nucleotides(sequence):
    """
    Find all dinucleotides and trinucleotides that exist in the sequence.
    Uses a single pass through the sequence (no brute force).

    Args:
        sequence: Input string sequence

    Returns:
        tuple: (set of dinucleotides, set of trinucleotides)
    """
    dinucleotides = set()
    trinucleotides = set()

   
    for i in range(len(sequence)):
        
        if i + 2 <= len(sequence):
            dinucleotides.add(sequence[i:i+2])

        
        if i + 3 <= len(sequence):
            trinucleotides.add(sequence[i:i+3])

    return dinucleotides, trinucleotides



if __name__ == "__main__":
    S = "ABAA"

    di, tri = find_nucleotides(S)

    print(f"Sequence: {S}")
    print(f"Dinucleotides: {sorted(di)}")
    print(f"Trinucleotides: {sorted(tri)}")

    
    print("\n" + "="*50 + "\n")

    dna_sequence = "ATCGATCG"
    di_dna, tri_dna = find_nucleotides(dna_sequence)

    print(f"DNA Sequence: {dna_sequence}")
    print(f"Dinucleotides: {sorted(di_dna)}")
    print(f"Trinucleotides: {sorted(tri_dna)}")
