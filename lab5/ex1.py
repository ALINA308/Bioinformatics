"""
DNA Sequence Reconstruction from Random Samples

This program:
1. Uses a real DNA sequence (1000-3000 nucleotides)
2. Takes 2000 random samples (100-150 bases each)
3. Reconstructs the original sequence using overlap assembly
"""

import random


def get_dna_sequence():
    sequence = """
ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGAC
CCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTAC
CTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGAC
CTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTG
GCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGC
TCCCTCTACCAGCTGGAGAACTACTGCAACTAGACGCAGCCCGCAGGCAGCCCCACACCC
GCCGCCTCCTGCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGCACCTGTGGGGCCGG
GAGGCTGGAAGCTGTGGCAGGAGCAGCCTCCAGGAGGCCAGGACAGCACAGCTTGGGGAC
CCTGGAGCAGCCCTGCAGCCCCTCTCCTTGGGGGCTGGAGGACAGAAGCAGGAGGCTGCA
GGCTCCAGGGCCAAGGAAGGATGGCCCCCAGGCAGCACAGCTCCCCCAGGATCAGCTGCC
TCCCGGAGGTGTGCGAGGAGCTGCCTGCCTCCCCGGGCTGCCCTGCCCCACTGGCCACTG
GGTGTCTGGAGCCAGGCATCTGCACCCGGAGCCTGGAGCCCTGGACTGCCTGCCTCCAGC
CTGCCCCTGCTCCCCCTGGCCTGGGGACCAGCTCAGAGCCTGGCCAGCTCCAGGTGCCCC
AGCCCTGGGCCTGGGACTCCCCTGCTGCCCCATCTGCCCCTTCTGCCTGCCCCCAGCCCA
CCCCCTGATGCCCTGGCCCCTTCCTGCAGCTCCCCCTCGGGGATGCTGCCCCAAATGTGG
GGACCTGGCAGAAGCTGGATCCCCAGCGCCCCCCAGCTGTGGCTGCAGAAGCTCCTCTTC
CCAGCTCCTGCCTCCTGCCCCTCCCTGCCCAGGCTCCTGCCCACCACCTTGGCCCTTGCT
GCCTTGGTTGCTGCCTTGGCCTTGGGCCAGGACTGGGGCCTGACCCCTGCCCTTTGCCCC
CAGACCTGTCCCACGTCTGTGGCTTTGCAGCCTGTGCCCTGCTCCCGCCCATGCTGGGTG
ATGCCGACACTGAGCTGAAGCTCCCCAAGTGCTGCCTGTGTGAGCTGGCCCTGCCTGGTG
CCCACAGACTGGGAGCTGCAGCCTGGGAGCTGCAGCCCCCGAGCTTGGGGCTGCAGCTCA
CCTCCAGCTGCTTTGACCACGCCTTGGCCCTGAGCTCCCCAGCACCTGGGCCAGCTGGGG
ACCTCCCAGGCCCAGAGCCCCGAGGAGCTGCTGGACGTGGCTGCCTTGGTTGTGGACCTG
GGCCAGGCTGGGGCCAGGGCCCAGCCCCCTGAGGGCCTGCTGGTGGTTGTGAGCCTGGGC
CAGGCTCCCAGCCCTGCTGGGACCTCGGGTTCCTGTCCGAGGCTGCCTTGGTGGTTGTGG
ACCTGGGCCAGGCTGGGGCTGAGGCTGGGGCCAGGTCCCCGGAGGGCTCCCCAGCTGCCC
TGGACCTGGGGGCTGAGGTCCTGAGCCCTGCTGGGGTGCTGAGGCTGCCCCTGGTTGTTG
TGGACCTGGGACAGGCTGGGACTGAGGCCAGCCCTGAGGGCCTGCTGGTGGTGGTGGACC
TGGGCCAGGCTCCCAGCTGCTGGGACCTCGGGTTCCTGTCCCAGGCTGCCTTGGTGGTTG
TGGACCTGGGCCAGGCTCCCAGCCCTGCTGGGACCTCGGGTTCCTGTCCGAGGCTGCCTT
GGTGGTTGTGGACCTGGGCCAGGCTGGGGCTGAGGCCGGGGCCAGGTCCCCGGAGGGCTC
CCCAGCAGCCCTGGACCTGGGGGCTGAGGTCCTGAGCCCTGCCGGGGTGCTGAGGCCGCC
CCTGGTTGTTGTGGACCTGGGACAGGCTGGGACTGAGGCTGGCCCTGAGGACCTGCCGGT
GGTGGTGGACCTG
""".replace('\n', '').replace(' ', '')

    return sequence[:2000]


def take_random_samples(sequence, num_samples=2000, min_length=100, max_length=150):
    samples = []
    seq_length = len(sequence)

    for _ in range(num_samples):
        sample_length = random.randint(min_length, max_length)

        if seq_length <= sample_length:
            start = 0
            sample_length = seq_length
        else:
            start = random.randint(0, seq_length - sample_length)

        sample = sequence[start:start + sample_length]
        samples.append(sample)

    return samples


def find_overlap(s1, s2, min_overlap=20):
    max_overlap = min(len(s1), len(s2))

    for overlap_len in range(max_overlap, min_overlap - 1, -1):
        if s1[-overlap_len:] == s2[:overlap_len]:
            return overlap_len

    return 0


def reconstruct_sequence(samples, min_overlap=20):
    if not samples:
        return ""

    unique_samples = list(set(samples))
    print(f"Original samples: {len(samples)}, Unique samples: {len(unique_samples)}")

    unique_samples.sort(key=len, reverse=True)
    reconstructed = unique_samples[0]
    remaining = unique_samples[1:]

    max_iterations = len(remaining) * 2
    iteration = 0

    while remaining and iteration < max_iterations:
        iteration += 1
        best_overlap = 0
        best_idx = -1
        best_position = None

        for i, sample in enumerate(remaining):
            overlap_end = find_overlap(reconstructed, sample, min_overlap)
            if overlap_end > best_overlap:
                best_overlap = overlap_end
                best_idx = i
                best_position = 'end'

            overlap_start = find_overlap(sample, reconstructed, min_overlap)
            if overlap_start > best_overlap:
                best_overlap = overlap_start
                best_idx = i
                best_position = 'start'

        if best_idx != -1:
            sample = remaining.pop(best_idx)
            if best_position == 'end':
                reconstructed += sample[best_overlap:]
            else:
                reconstructed = sample + reconstructed[best_overlap:]
        else:
            break

    print(f"Assembled using {len(unique_samples) - len(remaining)}/{len(unique_samples)} unique fragments")

    return reconstructed


def calculate_accuracy(original, reconstructed):
    stats = {
        'original_length': len(original),
        'reconstructed_length': len(reconstructed),
        'length_ratio': len(reconstructed) / len(original) if original else 0
    }

    max_match = 0

    for i in range(len(original)):
        for j in range(len(reconstructed)):
            k = 0
            while (i + k < len(original) and j + k < len(reconstructed) and
                   original[i + k] == reconstructed[j + k]):
                k += 1
            max_match = max(max_match, k)

    stats['longest_match'] = max_match
    stats['match_percentage'] = (max_match / len(original) * 100) if original else 0
    stats['exact_match'] = (original == reconstructed)

    return stats


def main():
    print("=" * 70)
    print("DNA SEQUENCE RECONSTRUCTION FROM RANDOM SAMPLES")
    print("=" * 70)

    print("\n[Step 1] Loading DNA sequence...")
    original_sequence = get_dna_sequence()
    print(f"Original sequence length: {len(original_sequence)} bp")
    print(f"First 100 bp: {original_sequence[:100]}...")

    print("\n[Step 2-3] Taking 2000 random samples (100-150 bp each)...")
    samples = take_random_samples(original_sequence, num_samples=2000,
                                  min_length=100, max_length=150)
    print(f"Total samples collected: {len(samples)}")
    print(f"Sample length range: {min(len(s) for s in samples)}-{max(len(s) for s in samples)} bp")

    print("\n[Step 4] Reconstructing sequence using overlap assembly...")
    reconstructed_sequence = reconstruct_sequence(samples, min_overlap=20)
    print(f"Reconstructed sequence length: {len(reconstructed_sequence)} bp")

    print("\n" + "=" * 70)
    print("RECONSTRUCTION RESULTS")
    print("=" * 70)

    stats = calculate_accuracy(original_sequence, reconstructed_sequence)

    print(f"Original length:        {stats['original_length']} bp")
    print(f"Reconstructed length:   {stats['reconstructed_length']} bp")
    print(f"Length ratio:           {stats['length_ratio']:.2%}")
    print(f"Longest match:          {stats['longest_match']} bp")
    print(f"Match percentage:       {stats['match_percentage']:.2f}%")
    print(f"Exact match:            {stats['exact_match']}")

    if stats['exact_match']:
        print("\n✓ SUCCESS: Perfectly reconstructed the original sequence!")
    else:
        print("\n✗ PARTIAL: Reconstruction differs from original")

    # Main problem analysis
    print("\n" + "=" * 70)
    print("MAIN PROBLEMS WITH THIS ALGORITHM APPROACH:")
    print("=" * 70)
    print("""
1. REPEAT SEQUENCES (Biggest Problem):
   - DNA sequences often contain repetitive regions
   - The greedy algorithm cannot distinguish which copy of a repeat to use
   - This leads to misassembly, collapsed repeats, or wrong connections
   - Example: ATGATGATGATG - which AT connects where?

2. AMBIGUOUS OVERLAPS:
   - Multiple fragments may have similar overlaps
   - Greedy algorithm picks first/best but may not be globally optimal
   - No backtracking if wrong choice is made

3. COVERAGE GAPS:
   - Random sampling may miss some regions entirely
   - Even with 2000 samples, some areas might have low/no coverage
   - Reconstruction will be incomplete

4. COMPUTATIONAL COMPLEXITY:
   - Greedy approach is O(n²) or worse for n fragments
   - Not scalable for large genomes
   - Modern assemblers use De Bruijn graphs (more efficient)

5. NO ERROR CORRECTION:
   - Doesn't handle sequencing errors
   - Single base difference breaks overlap detection
   - Real sequencing data has ~1% error rate

6. ORIENTATION:
   - This algorithm only checks forward direction
   - Real DNA fragments can be reverse-complemented
   - Would need to check both orientations

BETTER APPROACHES:
   - De Bruijn graph assembly (used in real assemblers)
   - Overlap-Layout-Consensus (OLC) with graph optimization
   - Use paired-end reads for scaffolding
   - Multiple sequence alignment
   - Machine learning-based methods
    """)


if __name__ == "__main__":
    main()
