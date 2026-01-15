"""
DNA Sequence Transition Matrix Calculator
Computes transition probabilities between nucleotides in a DNA sequence.
"""

import json
import os
import sys
from collections import defaultdict


class DNATransitionCalculator:
    """
    Calculates transition probabilities between DNA nucleotides.
    """

    def __init__(self, dna_sequence):
        """
        Initialize with a DNA sequence.

        Args:
            dna_sequence: String containing DNA sequence (A, C, G, T)
        """
        self.dna_sequence = dna_sequence.upper().strip()
        self.nucleotides = ['A', 'C', 'G', 'T']
        self.transition_counts = defaultdict(lambda: defaultdict(int))
        self.transition_matrix = None

        # Validate sequence
        self._validate_sequence()

    def _validate_sequence(self):
        """Validate that sequence contains only valid nucleotides."""
        valid_nucleotides = set(self.nucleotides)
        for nucleotide in self.dna_sequence:
            if nucleotide not in valid_nucleotides:
                raise ValueError(
                    f"Invalid nucleotide '{nucleotide}'. "
                    f"DNA sequence must contain only A, C, G, T"
                )

        if len(self.dna_sequence) < 2:
            raise ValueError("DNA sequence must be at least 2 nucleotides long")

    def calculate_transitions(self):
        """
        Count transitions between consecutive nucleotides.
        """
        # Count all transitions
        for i in range(len(self.dna_sequence) - 1):
            current = self.dna_sequence[i]
            next_nucleotide = self.dna_sequence[i + 1]
            self.transition_counts[current][next_nucleotide] += 1

        # Calculate probabilities
        self.transition_matrix = []

        for from_nuc in self.nucleotides:
            row = []
            total_transitions = sum(self.transition_counts[from_nuc].values())

            for to_nuc in self.nucleotides:
                if total_transitions > 0:
                    probability = self.transition_counts[from_nuc][to_nuc] / total_transitions
                else:
                    probability = 0.0
                row.append(round(probability, 4))

            self.transition_matrix.append(row)

        return self.transition_matrix

    def get_transition_counts(self):
        """Return the raw transition counts."""
        counts_dict = {}
        for from_nuc in self.nucleotides:
            counts_dict[from_nuc] = dict(self.transition_counts[from_nuc])
        return counts_dict

    def display_results(self):
        """Display the transition matrix and statistics."""
        print("=" * 70)
        print("DNA SEQUENCE TRANSITION MATRIX CALCULATOR")
        print("=" * 70)
        print(f"\nDNA Sequence ({len(self.dna_sequence)} nucleotides):")
        print(self.dna_sequence)
        print("\n" + "=" * 70)

        print("\nTransition Counts:")
        print("-" * 70)
        print(f"{'From/To':<10}", end="")
        for nuc in self.nucleotides:
            print(f"{nuc:>10}", end="")
        print()
        print("-" * 70)

        for from_nuc in self.nucleotides:
            print(f"{from_nuc:<10}", end="")
            for to_nuc in self.nucleotides:
                count = self.transition_counts[from_nuc][to_nuc]
                print(f"{count:>10}", end="")
            print()

        print("\n" + "=" * 70)
        print("\nTransition Probability Matrix:")
        print("-" * 70)
        print(f"{'From/To':<10}", end="")
        for nuc in self.nucleotides:
            print(f"{nuc:>10}", end="")
        print()
        print("-" * 70)

        for i, from_nuc in enumerate(self.nucleotides):
            print(f"{from_nuc:<10}", end="")
            for j, to_nuc in enumerate(self.nucleotides):
                prob = self.transition_matrix[i][j]
                print(f"{prob:>10.4f}", end="")
            print()

        print("\n" + "=" * 70)

    def save_to_json(self, output_file):
        """
        Save transition matrix and metadata to JSON file.

        Args:
            output_file: Path to output JSON file
        """
        output_data = {
            "dna_sequence": self.dna_sequence,
            "sequence_length": len(self.dna_sequence),
            "nucleotides": self.nucleotides,
            "transition_counts": self.get_transition_counts(),
            "transition_matrix": self.transition_matrix,
            "matrix_description": "Rows and columns ordered as: A, C, G, T"
        }

        with open(output_file, 'w') as f:
            json.dump(output_data, f, indent=2)

        print(f"Transition matrix saved to {output_file}")


def load_sequence_from_file(filename):
    """Load DNA sequence from a text file."""
    with open(filename, 'r') as f:
        sequence = f.read().strip()
    return sequence


def generate_random_dna_sequence(length=50):
    """Generate a random DNA sequence for testing."""
    import random
    nucleotides = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(nucleotides) for _ in range(length))


def main():
    """Main application entry point."""
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Check if sequence file is provided
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        if not os.path.isabs(input_file):
            if os.path.exists(input_file):
                input_file = os.path.abspath(input_file)
            else:
                input_file = os.path.join(script_dir, input_file)

        print(f"Loading DNA sequence from {input_file}...")
        dna_sequence = load_sequence_from_file(input_file)
    else:
        # Generate random sequence if no input provided
        print("No input file provided. Generating random DNA sequence...")
        print("Usage: python dna_transition.py [sequence_file.txt]\n")
        dna_sequence = generate_random_dna_sequence(50)

    # Create calculator and compute transitions
    calculator = DNATransitionCalculator(dna_sequence)
    calculator.calculate_transitions()
    calculator.display_results()

    # Determine output file
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
        if not os.path.isabs(output_file):
            output_dir = os.path.dirname(output_file)
            if output_dir and os.path.exists(output_dir):
                output_file = os.path.abspath(output_file)
            else:
                output_file = os.path.join(script_dir, output_file)
    else:
        output_file = os.path.join(script_dir, 'transition_matrix.json')

    # Save results
    calculator.save_to_json(output_file)


if __name__ == "__main__":
    main()
