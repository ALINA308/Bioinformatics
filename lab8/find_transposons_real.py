import json
from collections import defaultdict

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def find_inverted_repeats(sequence, min_length=4, max_length=6, min_spacing=10, max_spacing=100):
    inverted_repeats = []
    seq_upper = sequence.upper()

    for length in range(max_length, min_length - 1, -1):
        print(f"  Searching for {length} bp repeats...")
        for i in range(len(seq_upper) - length):
            left_seq = seq_upper[i:i+length]

            if 'N' in left_seq:
                continue

            rc = reverse_complement(left_seq)

            for j in range(i + length + min_spacing, min(i + length + max_spacing, len(seq_upper) - length + 1)):
                right_seq = seq_upper[j:j+length]

                if right_seq == rc:
                    inverted_repeats.append({
                        'left_start': i,
                        'left_end': i + length - 1,
                        'right_start': j,
                        'right_end': j + length - 1,
                        'sequence': left_seq,
                        'length': length,
                        'spacing': j - (i + length)
                    })

    filtered = []
    used = set()
    inverted_repeats.sort(key=lambda x: (-x['length'], x['left_start']))

    for ir in inverted_repeats:
        pos_range = set(range(ir['left_start'], ir['left_end'] + 1))
        pos_range.update(range(ir['right_start'], ir['right_end'] + 1))

        if not pos_range.intersection(used):
            filtered.append(ir)
            used.update(pos_range)

    return filtered

def read_fasta(filename):
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

def analyze_genome(filename, genome_name):
    print(f"\n{'='*70}")
    print(f"Analyzing {genome_name}")
    print(f"{'='*70}")

    sequence = read_fasta(filename)
    genome_length = len(sequence)
    print(f"Genome length: {genome_length:,} bp")

    print("\nSearching for inverted repeats (4-6 bp)...")
    inverted_repeats = find_inverted_repeats(sequence)

    print(f"Found {len(inverted_repeats)} potential transposon sites")

    stats = defaultdict(int)
    for ir in inverted_repeats:
        stats[ir['length']] += 1

    print("\nBreakdown by repeat length:")
    for length in sorted(stats.keys(), reverse=True):
        print(f"  {length} bp: {stats[length]} occurrences")

    top_repeats = inverted_repeats[:20]

    if top_repeats:
        print(f"\nTop 20 inverted repeats:")
        print(f"{'Pos':<12} {'Seq':<10} {'Len':<6} {'Spacing':<10} {'Region':<20}")
        print(f"{'-'*70}")

        for ir in top_repeats:
            region = f"{ir['left_start']}-{ir['right_end']}"
            print(f"{ir['left_start']:<12} {ir['sequence']:<10} {ir['length']:<6} {ir['spacing']:<10} {region:<20}")

    return {
        'genome_name': genome_name,
        'genome_length': genome_length,
        'total_inverted_repeats': len(inverted_repeats),
        'stats': dict(stats),
        'top_20': top_repeats
    }

def main():
    print("="*70)
    print("TRANSPOSON DETECTION IN BACTERIAL GENOMES")
    print("="*70)
    print("Searching for inverted repeats (4-6 bp) as markers of transposons")

    genomes = [
        ("Escherichia_coli.fasta", "Escherichia coli"),
        ("Bacillus_subtilis.fasta", "Bacillus subtilis"),
        ("Mycoplasma_genitalium.fasta", "Mycoplasma genitalium")
    ]

    results = []

    for filename, name in genomes:
        try:
            result = analyze_genome(filename, name)
            results.append(result)
        except FileNotFoundError:
            print(f"\nError: {filename} not found. Run download_genomes.py first.")

    if results:
        with open('bacterial_transposon_analysis.json', 'w') as f:
            json.dump(results, f, indent=2)

        print(f"\n{'='*70}")
        print("SUMMARY")
        print(f"{'='*70}")
        for r in results:
            print(f"\n{r['genome_name']}:")
            print(f"  Genome size: {r['genome_length']:,} bp")
            print(f"  Inverted repeats found: {r['total_inverted_repeats']}")
            density = (r['total_inverted_repeats'] / r['genome_length']) * 1000000
            print(f"  Density: {density:.2f} per Mbp")

        print(f"\nResults saved to bacterial_transposon_analysis.json")
        print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
