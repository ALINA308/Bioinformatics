"""
DNA Repetition Detector
Detects repetitions of 6-10 base pairs in a DNA sequence
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_dna_sequence(filename):
    sequence = ""
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('>'):
                sequence += line.upper()
    return sequence


def detect_repetitions(sequence, min_length=6, max_length=10):
    repetitions = {}

    for pattern_length in range(min_length, max_length + 1):
        for i in range(len(sequence) - pattern_length + 1):
            pattern = sequence[i:i + pattern_length]

            positions = []
            for j in range(len(sequence) - pattern_length + 1):
                if sequence[j:j + pattern_length] == pattern:
                    positions.append(j)

            if len(positions) > 1:
                if pattern not in repetitions:
                    repetitions[pattern] = positions

    return repetitions


def filter_repetitions(repetitions, min_occurrences=2):
    return {pattern: positions for pattern, positions in repetitions.items()
            if len(positions) >= min_occurrences}


def display_results(repetitions, sequence, top_n=50):
    print(f"\n{'='*80}")
    print(f"DNA REPETITION ANALYSIS REPORT")
    print(f"{'='*80}")
    print(f"\nSequence length: {len(sequence)} base pairs")
    print(f"Total unique repetitive patterns found: {len(repetitions)}")

    sorted_repetitions = sorted(repetitions.items(),
                               key=lambda x: len(x[1]),
                               reverse=True)

    print(f"\n{'='*80}")
    print(f"TOP {min(top_n, len(sorted_repetitions))} MOST FREQUENT REPETITIONS")
    print(f"{'='*80}")
    print(f"\n{'Pattern':<15} {'Length':<8} {'Count':<8} {'Positions'}")
    print(f"{'-'*80}")

    for i, (pattern, positions) in enumerate(sorted_repetitions[:top_n]):
        if len(positions) > 10:
            pos_str = ', '.join(map(str, positions[:10])) + f'... (+{len(positions)-10} more)'
        else:
            pos_str = ', '.join(map(str, positions))

        print(f"{pattern:<15} {len(pattern):<8} {len(positions):<8} {pos_str}")

    print(f"\n{'='*80}")
    print(f"STATISTICS BY PATTERN LENGTH")
    print(f"{'='*80}")
    print(f"\n{'Length':<10} {'Unique Patterns':<20} {'Total Occurrences'}")
    print(f"{'-'*80}")

    for length in range(6, 11):
        patterns_of_length = {p: pos for p, pos in repetitions.items() if len(p) == length}
        total_occurrences = sum(len(pos) for pos in patterns_of_length.values())
        print(f"{length:<10} {len(patterns_of_length):<20} {total_occurrences}")


def plot_repetition_frequencies(repetitions):
    sorted_repetitions = sorted(repetitions.items(), key=lambda x: len(x[1]), reverse=True)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

    # Plot 1: Top 20 most frequent patterns
    top_20 = sorted_repetitions[:20]
    patterns = [p[0] for p in top_20]
    frequencies = [len(p[1]) for p in top_20]

    ax1.barh(range(len(patterns)), frequencies, color='steelblue')
    ax1.set_yticks(range(len(patterns)))
    ax1.set_yticklabels(patterns, fontsize=8)
    ax1.set_xlabel('Frequency (Number of Occurrences)', fontsize=10)
    ax1.set_title('Top 20 Most Frequent Patterns', fontsize=12, fontweight='bold')
    ax1.invert_yaxis()
    ax1.grid(axis='x', alpha=0.3)

    # Plot 2: Frequency distribution histogram
    all_frequencies = [len(positions) for positions in repetitions.values()]
    ax2.hist(all_frequencies, bins=30, color='coral', edgecolor='black', alpha=0.7)
    ax2.set_xlabel('Number of Occurrences', fontsize=10)
    ax2.set_ylabel('Number of Patterns', fontsize=10)
    ax2.set_title('Frequency Distribution of Patterns', fontsize=12, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3)

    # Plot 3: Statistics by pattern length
    lengths = list(range(6, 11))
    unique_patterns = []
    total_occurrences = []

    for length in lengths:
        patterns_of_length = {p: pos for p, pos in repetitions.items() if len(p) == length}
        unique_patterns.append(len(patterns_of_length))
        total_occurrences.append(sum(len(pos) for pos in patterns_of_length.values()))

    x = range(len(lengths))
    width = 0.35

    ax3.bar([i - width/2 for i in x], unique_patterns, width, label='Unique Patterns', color='mediumseagreen')
    ax3.bar([i + width/2 for i in x], total_occurrences, width, label='Total Occurrences', color='mediumpurple')
    ax3.set_xlabel('Pattern Length (base pairs)', fontsize=10)
    ax3.set_ylabel('Count', fontsize=10)
    ax3.set_title('Statistics by Pattern Length', fontsize=12, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels(lengths)
    ax3.legend()
    ax3.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig('repetition_frequency_plots.png', dpi=300, bbox_inches='tight')
    print("Plots saved to 'repetition_frequency_plots.png'")
    plt.close()


def save_results_to_file(repetitions, sequence, output_file):
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("DNA REPETITION ANALYSIS REPORT\n")
        f.write("="*80 + "\n\n")
        f.write(f"Sequence length: {len(sequence)} base pairs\n")
        f.write(f"Total unique repetitive patterns found: {len(repetitions)}\n\n")

        sorted_repetitions = sorted(repetitions.items(),
                                   key=lambda x: len(x[1]),
                                   reverse=True)

        f.write("="*80 + "\n")
        f.write("ALL REPETITIVE PATTERNS (sorted by frequency)\n")
        f.write("="*80 + "\n\n")
        f.write(f"{'Pattern':<15} {'Length':<8} {'Count':<8} {'Positions'}\n")
        f.write("-"*80 + "\n")

        for pattern, positions in sorted_repetitions:
            pos_str = ', '.join(map(str, positions))
            f.write(f"{pattern:<15} {len(pattern):<8} {len(positions):<8} {pos_str}\n")

        f.write("\n" + "="*80 + "\n")
        f.write("STATISTICS BY PATTERN LENGTH\n")
        f.write("="*80 + "\n\n")
        f.write(f"{'Length':<10} {'Unique Patterns':<20} {'Total Occurrences'}\n")
        f.write("-"*80 + "\n")

        for length in range(6, 11):
            patterns_of_length = {p: pos for p, pos in repetitions.items() if len(p) == length}
            total_occurrences = sum(len(pos) for pos in patterns_of_length.values())
            f.write(f"{length:<10} {len(patterns_of_length):<20} {total_occurrences}\n")


def main():
    input_file = "dna_sequence.txt"
    output_file = "repetition_analysis.txt"

    print("Loading DNA sequence...")
    sequence = read_dna_sequence(input_file)

    print(f"Sequence loaded: {len(sequence)} base pairs")
    print("\nDetecting repetitions (6-10 base pairs)...")

    repetitions = detect_repetitions(sequence, min_length=6, max_length=10)
    repetitions = filter_repetitions(repetitions, min_occurrences=2)

    display_results(repetitions, sequence)

    print(f"\n{'='*80}")
    print(f"Saving detailed results to '{output_file}'...")
    save_results_to_file(repetitions, sequence, output_file)
    print(f"Results saved successfully!")
    print(f"{'='*80}\n")

    print("Generating frequency plots...")
    plot_repetition_frequencies(repetitions)


if __name__ == "__main__":
    main()
