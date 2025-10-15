"""
Sliding Window Melting Temperature Calculator
Calculates melting temperature (Tm) over a DNA sequence using a sliding window approach.
Window size: 8 base pairs
Input: FASTA file
Output: Chart/plot of melting temperatures
"""

import matplotlib.pyplot as plt


def read_fasta(filename):
    
    header = ""
    sequence = ""

    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    header = line[1:]  # Remove '>' character
                else:
                    sequence += line.upper()

        if not sequence:
            raise ValueError("No sequence found in FASTA file")

        return header, sequence

    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file '{filename}' not found")
    except Exception as e:
        raise Exception(f"Error reading FASTA file: {str(e)}")


def calculate_melting_temperature(sequence):
    
    sequence = sequence.upper()

    # Count nucleotides
    count_A = sequence.count('A')
    count_T = sequence.count('T')
    count_G = sequence.count('G')
    count_C = sequence.count('C')

    length = len(sequence)

    if length == 0:
        return 0.0

    if length <= 14:
        tm = 2 * (count_A + count_T) + 4 * (count_G + count_C)
    else:
        gc_content = count_G + count_C
        tm = 64.9 + 41 * (gc_content - 16.4) / length

    return tm


def sliding_window_tm(sequence, window_size=8):
   
    results = []
    sequence = sequence.upper()
    seq_length = len(sequence)

    if seq_length < window_size:
        print(f"Warning: Sequence length ({seq_length}) is shorter than window size ({window_size})")
        print(f"Calculating Tm for the entire sequence")
        tm = calculate_melting_temperature(sequence)
        results.append((0, sequence, tm))
        return results

    for i in range(seq_length - window_size + 1):
        window = sequence[i:i + window_size]
        tm = calculate_melting_temperature(window)
        results.append((i, window, tm))

    return results


def display_results(header, sequence, results, window_size):
   
    print("\n" + "="*70)
    print("MELTING TEMPERATURE ANALYSIS - SLIDING WINDOW METHOD")
    print("="*70)
    print(f"\nSequence Header: {header}")
    print(f"Sequence Length: {len(sequence)} bp")
    print(f"Window Size: {window_size} bp")
    print(f"Number of Windows: {len(results)}")
    print("\n" + "-"*70)
    print(f"{'Position':<10} {'Window Sequence':<15} {'Tm (C)':<10}")
    print("-"*70)

    for pos, window, tm in results:
        print(f"{pos:<10} {window:<15} {tm:>6.2f}")

    if results:
        tm_values = [tm for _, _, tm in results]
        avg_tm = sum(tm_values) / len(tm_values)
        min_tm = min(tm_values)
        max_tm = max(tm_values)

        print("\n" + "-"*70)
        print("STATISTICS:")
        print(f"  Average Tm: {avg_tm:.2f} C")
        print(f"  Minimum Tm: {min_tm:.2f} C")
        print(f"  Maximum Tm: {max_tm:.2f} C")
        print("="*70 + "\n")


def plot_melting_temperatures(header, results, window_size):
    
    if not results:
        print("No data to plot.")
        return

    
    positions = [pos for pos, _, _ in results]
    tm_values = [tm for _, _, tm in results]

  
    avg_tm = sum(tm_values) / len(tm_values)

   
    plt.figure(figsize=(12, 6))
    plt.plot(positions, tm_values, 'b-', linewidth=2, label='Melting Temperature')
    plt.axhline(y=avg_tm, color='r', linestyle='--', linewidth=1, label=f'Average Tm: {avg_tm:.2f}°C')

   
    plt.xlabel('Position in Sequence (bp)', fontsize=12)
    plt.ylabel('Melting Temperature (°C)', fontsize=12)
    plt.title(f'Melting Temperature Profile - Sliding Window (size: {window_size} bp)\n{header}', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()

   
    y_margin = (max(tm_values) - min(tm_values)) * 0.1
    plt.ylim(min(tm_values) - y_margin, max(tm_values) + y_margin)

    plt.tight_layout()
    plt.show()

    print("\nChart displayed successfully!")


def main():
   
    import sys

    print("\n" + "="*70)
    print("MELTING TEMPERATURE CALCULATOR - SLIDING WINDOW METHOD")
    print("="*70)

   
    if len(sys.argv) > 1:
        fasta_file = sys.argv[1]
    else:
        fasta_file = input("\nEnter the path to the FASTA file: ").strip()

    try:
       
        print(f"\nReading FASTA file: {fasta_file}")
        header, sequence = read_fasta(fasta_file)

       
        window_size = 8

        print(f"Calculating melting temperatures with window size {window_size}...")
        results = sliding_window_tm(sequence, window_size)

        display_results(header, sequence, results, window_size)

       
        print("\nGenerating chart...")
        plot_melting_temperatures(header, results, window_size)

    except Exception as e:
        print(f"\nError: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
