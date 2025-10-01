import sys
import os
from collections import Counter

# Problem 1: Find alphabet of a sequence
def find_alphabet(sequence: str) -> set:
    return set(sequence)

# Problem 2: DNA sequence analysis
seq = "ACGGGCATATGCGC"
alphabet = find_alphabet(seq)
print("Problem 1 - Alphabet of sequence:", alphabet)

print("\nProblem 2 - DNA Sequence Analysis:")
print(f"Sequence: {seq}")
percentages = {base: (seq.count(base) / len(seq)) * 100 for base in alphabet}
for base in sorted(percentages):
    print(f"{base}: {percentages[base]:.2f}%")

# Problem 3: FASTA file reader
def read_fasta_file(file_path):
    
    sequences = []
    current_sequence = ""
    current_header = ""

    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if current_sequence:
                        sequences.append((current_header, current_sequence))
                    current_header = line[1:]
                    current_sequence = ""
                else:
                    current_sequence += line.upper()

            if current_sequence:
                sequences.append((current_header, current_sequence))

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

    return sequences

def analyze_sequence_composition(sequence):
    
    if not sequence:
        return {}

    alphabet = find_alphabet(sequence)
    percentages = {base: (sequence.count(base) / len(sequence)) * 100 for base in alphabet}
    return percentages

def display_fasta_results(sequences_data):
    
    for i, (header, sequence) in enumerate(sequences_data, 1):
        print(f"\n{'='*60}")
        print(f"Sequence {i}: {header}")
        print(f"{'='*60}")
        print(f"Length: {len(sequence)} characters")
        print(f"Alphabet: {find_alphabet(sequence)}")

        composition = analyze_sequence_composition(sequence)

        if composition:
            print("\nComposition Analysis:")
            print("-" * 40)
            for base in sorted(composition):
                print(f"{base}: {composition[base]:.2f}%")

        print(f"\nSequence preview: {sequence[:50]}{'...' if len(sequence) > 50 else ''}")

print("\n" + "="*60)
print("Problem 3 - FASTA File Analysis")
print("="*60)

if len(sys.argv) > 1:
    file_path = sys.argv[1]
    if os.path.exists(file_path):
        print(f"\nReading FASTA file: {file_path}")
        sequences = read_fasta_file(file_path)

        if sequences and len(sequences) > 0:
            print(f"\nFound {len(sequences)} sequence(s)")
            display_fasta_results(sequences)
        else:
            print("No sequences found in the file.")
    else:
        print(f"Error: File '{file_path}' does not exist.")
else:
    print("To analyze a FASTA file, run: python ex1_2py filename.fasta")
    

    

