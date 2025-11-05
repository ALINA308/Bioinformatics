import random
import matplotlib.pyplot as plt
import os

def read_fasta_sequence(filename):
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

def extract_random_samples(sequence, num_samples=10, min_length=100, max_length=3000):
    samples = []
    seq_length = len(sequence)

    for i in range(num_samples):
        fragment_length = random.randint(min_length, min(max_length, seq_length))

        if seq_length > fragment_length:
            start_pos = random.randint(0, seq_length - fragment_length)
        else:
            start_pos = 0
            fragment_length = seq_length

        fragment = sequence[start_pos:start_pos + fragment_length]

        samples.append({
            'id': i + 1,
            'sequence': fragment,
            'length': len(fragment),
            'start_pos': start_pos
        })

    return samples

def simulate_gel_electrophoresis(samples):
    max_length = max(sample['length'] for sample in samples)

    for sample in samples:
        migration_factor = 1 - (sample['length'] / max_length)
        sample['migration_distance'] = migration_factor

    return samples

def visualize_gel(samples, ladder_sizes=[500, 1500, 3000]):
    fig, ax = plt.subplots(figsize=(14, 8))

    ax.set_facecolor('black')
    fig.patch.set_facecolor('white')

    max_bp = max(sample['length'] for sample in samples)
    if max_bp > max(ladder_sizes):
        max_ladder = max_bp * 1.05
        if max_ladder not in ladder_sizes:
            ladder_sizes = sorted(ladder_sizes + [int(max_ladder)])
    else:
        max_ladder = max(ladder_sizes)

    lane_width = 0.7
    lane_x = 0

    print(f"\nLadder scale (max = {max_ladder:.0f} bp):")
    for size in ladder_sizes:
        y_pos = size / max_ladder
        print(f"  {size} bp -> y_pos = {y_pos:.3f}")
        ax.barh(y_pos, lane_width, height=0.015, left=lane_x - lane_width/2,
                color='white', edgecolor='gray', linewidth=0.8)

    left_label_x = -0.5
    for size in ladder_sizes:
        y_pos = size / max_ladder
        ax.text(left_label_x, y_pos, f'{size} bp',
                ha='right', va='center', fontsize=14, color='red', fontweight='bold')

    print("\nSample positioning (bp value -> position):")
    for idx, sample in enumerate(samples, start=1):
        lane_x = idx * 1.1
        y_pos = sample['length'] / max_ladder
        print(f"  Sample {idx}: {sample['length']} bp -> y_pos = {y_pos:.3f} (out of 1.0)")

        ax.barh(y_pos, lane_width, height=0.015, left=lane_x - lane_width/2,
                color='white', edgecolor='gray', linewidth=0.8, alpha=1.0)

    right_edge = (len(samples)) * 1.1 + 1.2
    for sample in samples:
        y_pos = sample['length'] / max_ladder
        ax.text(right_edge, y_pos, f'{sample["length"]} bp',
                ha='left', va='center', fontsize=11, color='black', fontweight='bold')

    ax.set_xlim(-1.0, right_edge + 1.0)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel('Lanes', fontsize=12, fontweight='bold')
    ax.set_ylabel('Migration Distance →', fontsize=12, fontweight='bold')
    ax.set_title('Gel Electrophoresis Simulation\n(Shorter fragments migrate further)',
                 fontsize=14, fontweight='bold', pad=20)

    lane_labels = ['Ladder'] + [f'S{i}' for i in range(1, len(samples) + 1)]
    lane_positions = [0] + [i * 1.1 for i in range(1, len(samples) + 1)]
    ax.set_xticks(lane_positions)
    ax.set_xticklabels(lane_labels)

    ax.set_yticks([])

    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)

    for lane_x in lane_positions:
        well = plt.Rectangle((lane_x - lane_width/2, 0.95), lane_width, 0.05,
                            facecolor='#404040', edgecolor='white', linewidth=2)
        ax.add_patch(well)

    plt.tight_layout()
    plt.savefig('c:/Users/alina/Desktop/AN4/bioinformatics/lab6/gel_electrophoresis.png',
                dpi=300, bbox_inches='tight', pad_inches=0.3, facecolor='white')
    print("\n> Gel image saved as 'gel_electrophoresis.png'")
    plt.close()

def print_results(original_sequence, samples):
    print("\n" + "="*75)
    print("DNA FRAGMENT ANALYSIS RESULTS")
    print("="*75)
    print(f"Original sequence length: {len(original_sequence)} bp")
    print(f"Number of samples: {len(samples)}")
    print("\n" + "-"*75)
    print(f"{'Sample':<10} {'Length (bp)':<15} {'Start Pos':<15} {'Migration':<15}")
    print("-"*75)

    sorted_samples = sorted(samples, key=lambda x: x['migration_distance'], reverse=True)

    for sample in sorted_samples:
        print(f"Sample {sample['id']:<3} {sample['length']:<15} {sample['start_pos']:<15} {sample['migration_distance']:.4f}")

    print("="*75)
    print("\nNote: Higher migration distance = shorter fragment = travels further")

def main():
    print("="*75)
    print(" GEL ELECTROPHORESIS SIMULATION")
    print("="*75)

    print("\n[1] Reading DNA sequence from FASTA file...")
    sequence = read_fasta_sequence('c:/Users/alina/Desktop/AN4/bioinformatics/lab6/sequence.fasta')
    print(f"    > Sequence loaded: {len(sequence)} bp")
    print(f"    First 60 bp: {sequence[:60]}...")

    print("\n[2] Extracting 10 random samples from the sequence...")
    samples = extract_random_samples(sequence, num_samples=10, min_length=100, max_length=len(sequence))
    print(f"    > {len(samples)} samples extracted and stored in array")

    print("\n[3] Simulating gel electrophoresis migration...")
    samples = simulate_gel_electrophoresis(samples)
    print("    > Migration simulation complete")

    print_results(sequence, samples)

    print("\n[4] Creating visual representation of the gel...")
    visualize_gel(samples)

    print("\n[5] Opening the gel image...")
    image_path = 'c:/Users/alina/Desktop/AN4/bioinformatics/lab6/gel_electrophoresis.png'
    os.startfile(image_path)

    print("\n" + "="*75)
    print("SIMULATION COMPLETE!")
    print("="*75)

if __name__ == "__main__":
    main()
