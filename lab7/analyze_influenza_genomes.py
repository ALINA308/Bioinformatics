"""
Influenza Genomes Repetition Analysis
Analyzes 10 influenza virus genomes and plots repetition frequencies for each
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

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


def plot_genome_frequencies(repetitions, genome_name, output_file):
    sorted_repetitions = sorted(repetitions.items(), key=lambda x: len(x[1]), reverse=True)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'Repetition Analysis: {genome_name}', fontsize=14, fontweight='bold')

    # Plot 1: Top 20 most frequent patterns
    top_20 = sorted_repetitions[:20]
    if top_20:
        patterns = [p[0] for p in top_20]
        frequencies = [len(p[1]) for p in top_20]

        ax1.barh(range(len(patterns)), frequencies, color='steelblue')
        ax1.set_yticks(range(len(patterns)))
        ax1.set_yticklabels(patterns, fontsize=8)
        ax1.set_xlabel('Frequency (Number of Occurrences)', fontsize=10)
        ax1.set_title('Top 20 Most Frequent Patterns', fontsize=11, fontweight='bold')
        ax1.invert_yaxis()
        ax1.grid(axis='x', alpha=0.3)

    # Plot 2: Frequency distribution histogram
    all_frequencies = [len(positions) for positions in repetitions.values()]
    ax2.hist(all_frequencies, bins=20, color='coral', edgecolor='black', alpha=0.7)
    ax2.set_xlabel('Number of Occurrences', fontsize=10)
    ax2.set_ylabel('Number of Patterns', fontsize=10)
    ax2.set_title('Frequency Distribution', fontsize=11, fontweight='bold')
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
    ax3.set_xlabel('Pattern Length (bp)', fontsize=10)
    ax3.set_ylabel('Count', fontsize=10)
    ax3.set_title('Statistics by Pattern Length', fontsize=11, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels(lengths)
    ax3.legend()
    ax3.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def analyze_genome(genome_file, genome_name, output_dir):
    print(f"\nAnalyzing {genome_name}...")

    sequence = read_dna_sequence(genome_file)
    print(f"  Sequence length: {len(sequence)} bp")

    repetitions = detect_repetitions(sequence, min_length=6, max_length=10)
    repetitions = filter_repetitions(repetitions, min_occurrences=2)

    print(f"  Unique repetitive patterns found: {len(repetitions)}")

    if repetitions:
        sorted_reps = sorted(repetitions.items(), key=lambda x: len(x[1]), reverse=True)
        top_pattern, top_positions = sorted_reps[0]
        print(f"  Most frequent pattern: {top_pattern} ({len(top_positions)} occurrences)")

    plot_file = os.path.join(output_dir, f'{genome_name}_frequency_plot.png')
    plot_genome_frequencies(repetitions, genome_name, plot_file)
    print(f"  Plot saved: {plot_file}")

    return {
        'name': genome_name,
        'length': len(sequence),
        'total_patterns': len(repetitions),
        'repetitions': repetitions
    }


def create_summary_comparison(results, output_file):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Influenza Genomes Comparison - Repetition Analysis', fontsize=16, fontweight='bold')

    genome_names = [r['name'] for r in results]

    # Plot 1: Total unique patterns per genome
    total_patterns = [r['total_patterns'] for r in results]
    ax1.bar(range(len(genome_names)), total_patterns, color='steelblue', alpha=0.7)
    ax1.set_xlabel('Genome', fontsize=10)
    ax1.set_ylabel('Number of Unique Patterns', fontsize=10)
    ax1.set_title('Total Unique Repetitive Patterns', fontsize=12, fontweight='bold')
    ax1.set_xticks(range(len(genome_names)))
    ax1.set_xticklabels([f'G{i+1}' for i in range(len(genome_names))], fontsize=9)
    ax1.grid(axis='y', alpha=0.3)

    # Plot 2: Genome lengths
    lengths = [r['length'] for r in results]
    ax2.bar(range(len(genome_names)), lengths, color='coral', alpha=0.7)
    ax2.set_xlabel('Genome', fontsize=10)
    ax2.set_ylabel('Length (base pairs)', fontsize=10)
    ax2.set_title('Genome Sequence Lengths', fontsize=12, fontweight='bold')
    ax2.set_xticks(range(len(genome_names)))
    ax2.set_xticklabels([f'G{i+1}' for i in range(len(genome_names))], fontsize=9)
    ax2.grid(axis='y', alpha=0.3)

    # Plot 3: Pattern distribution by length
    pattern_lengths = {length: [] for length in range(6, 11)}
    for result in results:
        for length in range(6, 11):
            count = sum(1 for p in result['repetitions'].keys() if len(p) == length)
            pattern_lengths[length].append(count)

    x = range(len(genome_names))
    width = 0.15
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7']

    for i, (length, counts) in enumerate(pattern_lengths.items()):
        offset = width * (i - 2)
        ax3.bar([xi + offset for xi in x], counts, width, label=f'{length} bp', color=colors[i], alpha=0.8)

    ax3.set_xlabel('Genome', fontsize=10)
    ax3.set_ylabel('Number of Patterns', fontsize=10)
    ax3.set_title('Pattern Distribution by Length', fontsize=12, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels([f'G{i+1}' for i in range(len(genome_names))], fontsize=9)
    ax3.legend(title='Length', fontsize=8)
    ax3.grid(axis='y', alpha=0.3)

    # Plot 4: Top pattern frequency comparison
    max_frequencies = []
    for result in results:
        if result['repetitions']:
            max_freq = max(len(positions) for positions in result['repetitions'].values())
            max_frequencies.append(max_freq)
        else:
            max_frequencies.append(0)

    ax4.bar(range(len(genome_names)), max_frequencies, color='mediumseagreen', alpha=0.7)
    ax4.set_xlabel('Genome', fontsize=10)
    ax4.set_ylabel('Maximum Pattern Frequency', fontsize=10)
    ax4.set_title('Highest Pattern Repetition Count', fontsize=12, fontweight='bold')
    ax4.set_xticks(range(len(genome_names)))
    ax4.set_xticklabels([f'G{i+1}' for i in range(len(genome_names))], fontsize=9)
    ax4.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    genomes_dir = 'influenza_genomes'
    output_dir = 'influenza_analysis_results'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("="*80)
    print("INFLUENZA GENOMES REPETITION ANALYSIS")
    print("="*80)

    results = []

    for i in range(1, 11):
        genome_file = os.path.join(genomes_dir, f'genome_{i}.txt')
        genome_name = f'Genome_{i}'

        result = analyze_genome(genome_file, genome_name, output_dir)
        results.append(result)

    print("\n" + "="*80)
    print("Creating comparison summary...")
    summary_file = os.path.join(output_dir, 'genomes_comparison_summary.png')
    create_summary_comparison(results, summary_file)
    print(f"Summary comparison saved: {summary_file}")

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"\nTotal genomes analyzed: {len(results)}")
    print(f"\nIndividual genome statistics:")
    for i, result in enumerate(results, 1):
        print(f"  Genome {i}: {result['length']} bp, {result['total_patterns']} unique patterns")

    print("\n" + "="*80)
    print("Analysis complete! Check the 'influenza_analysis_results' folder for:")
    print("  - Individual genome frequency plots")
    print("  - Comparative summary plot")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
