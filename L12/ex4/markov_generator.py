"""
Markov Chain-Based Sequence Generator for Bioinformatics
Generates new sequences based on learned transition probabilities.
"""

import json
import os
import sys
import numpy as np
from collections import defaultdict, Counter


class MarkovSequenceGenerator:
    """
    A Markov chain-based sequence generator that uses transition probabilities
    to generate new sequences.
    """

    def __init__(self, transition_matrix):
        """
        Initialize the Markov chain generator.

        Args:
            transition_matrix: Dictionary mapping states to their transition probabilities
                              Format: {"A": {"A": 0.5, "B": 0.5}, "B": {...}, ...}
        """
        self.transition_matrix = transition_matrix
        self.states = list(transition_matrix.keys())

        # Validate the matrix
        self._validate_matrix()

        # Prepare data structures for efficient sampling
        self.transition_probs = {}
        self.transition_states = {}
        self._prepare_sampling_structures()

    def _validate_matrix(self):
        """
        Validate that the transition matrix is properly formed.
        - All probabilities are non-negative
        - Probabilities sum to 1.0 for each state (with tolerance)
        """
        if not self.states:
            raise ValueError("Transition matrix is empty")

        for state in self.states:
            transitions = self.transition_matrix[state]

            if not transitions:
                raise ValueError(f"State '{state}' has no transitions")

            # Check probabilities are non-negative
            for next_state, prob in transitions.items():
                if prob < 0:
                    raise ValueError(
                        f"Negative probability {prob} for transition {state}->{next_state}"
                    )

            # Check probabilities sum to 1.0 (with small tolerance for floating point)
            total_prob = sum(transitions.values())
            if abs(total_prob - 1.0) > 1e-6:
                raise ValueError(
                    f"Probabilities for state '{state}' sum to {total_prob}, not 1.0"
                )

    def _prepare_sampling_structures(self):
        """
        Pre-compute arrays for efficient random sampling.
        For each state, store the list of possible next states and their probabilities.
        """
        for state in self.states:
            transitions = self.transition_matrix[state]
            next_states = list(transitions.keys())
            probs = [transitions[s] for s in next_states]

            self.transition_states[state] = next_states
            self.transition_probs[state] = probs

    def generate_sequence(self, length, start_state=None):
        """
        Generate a sequence of specified length using the Markov chain.

        Args:
            length: Number of characters to generate
            start_state: Starting state (if None, choose randomly)

        Returns:
            Generated sequence as a string
        """
        if length <= 0:
            raise ValueError("Sequence length must be positive")

        # Choose starting state
        if start_state is None:
            current_state = np.random.choice(self.states)
        else:
            if start_state not in self.states:
                raise ValueError(f"Unknown starting state: {start_state}")
            current_state = start_state

        sequence = [current_state]

        # Generate remaining characters
        for _ in range(length - 1):
            # Get possible next states and their probabilities
            next_states = self.transition_states[current_state]
            probs = self.transition_probs[current_state]

            # Sample next state based on transition probabilities
            current_state = np.random.choice(next_states, p=probs)
            sequence.append(current_state)

        return ''.join(sequence)

    def generate_multiple_sequences(self, num_sequences, lengths):
        """
        Generate multiple sequences of varying lengths.

        Args:
            num_sequences: Number of sequences to generate
            lengths: List of lengths or single length for all sequences

        Returns:
            List of generated sequences
        """
        if isinstance(lengths, int):
            lengths = [lengths] * num_sequences

        if len(lengths) != num_sequences:
            raise ValueError(
                f"Number of lengths ({len(lengths)}) doesn't match num_sequences ({num_sequences})"
            )

        sequences = []
        for length in lengths:
            seq = self.generate_sequence(length)
            sequences.append(seq)

        return sequences

    def calculate_character_frequencies(self, sequences):
        """
        Calculate character frequency distribution in generated sequences.

        Args:
            sequences: List of sequences

        Returns:
            Dictionary of character frequencies
        """
        all_chars = ''.join(sequences)
        total_chars = len(all_chars)

        counter = Counter(all_chars)
        frequencies = {char: count / total_chars for char, count in counter.items()}

        return frequencies, counter

    def calculate_empirical_transitions(self, sequences):
        """
        Calculate empirical transition probabilities from generated sequences.

        Args:
            sequences: List of sequences

        Returns:
            Dictionary of empirical transition probabilities
        """
        transition_counts = defaultdict(lambda: defaultdict(int))
        state_counts = defaultdict(int)

        # Count transitions
        for sequence in sequences:
            for i in range(len(sequence) - 1):
                current = sequence[i]
                next_state = sequence[i + 1]
                transition_counts[current][next_state] += 1
                state_counts[current] += 1

        # Calculate probabilities
        empirical_matrix = {}
        for state in transition_counts:
            empirical_matrix[state] = {}
            total = state_counts[state]
            for next_state, count in transition_counts[state].items():
                empirical_matrix[state][next_state] = count / total

        return empirical_matrix

    def compare_matrices(self, empirical_matrix):
        """
        Compare empirical transition probabilities with the original matrix.

        Args:
            empirical_matrix: Empirical transition matrix from generated sequences

        Returns:
            Dictionary with comparison statistics
        """
        differences = {}

        for state in self.states:
            if state not in empirical_matrix:
                continue

            differences[state] = {}

            # Get all possible next states from both matrices
            all_next_states = set(self.transition_matrix[state].keys()) | set(empirical_matrix[state].keys())

            for next_state in all_next_states:
                expected = self.transition_matrix[state].get(next_state, 0.0)
                observed = empirical_matrix[state].get(next_state, 0.0)
                diff = abs(expected - observed)
                differences[state][next_state] = {
                    'expected': expected,
                    'observed': observed,
                    'difference': diff
                }

        return differences

    def display_results(self, sequences, show_full_comparison=False):
        """
        Display comprehensive results of sequence generation.

        Args:
            sequences: List of generated sequences
            show_full_comparison: Whether to show full transition comparison
        """
        print("=" * 80)
        print("MARKOV CHAIN SEQUENCE GENERATOR - RESULTS")
        print("=" * 80)

        # Display original transition matrix
        print("\nOriginal Transition Matrix:")
        print("-" * 80)
        for state in sorted(self.states):
            print(f"\nFrom '{state}':")
            for next_state, prob in sorted(self.transition_matrix[state].items()):
                print(f"  -> '{next_state}': {prob:.4f}")

        print("\n" + "=" * 80)
        print("\nGenerated Sequences:")
        print("-" * 80)
        for i, seq in enumerate(sequences, 1):
            print(f"{i:2d}. (length={len(seq):3d}) {seq}")

        # Calculate character frequencies
        print("\n" + "=" * 80)
        print("\nCharacter Frequency Statistics:")
        print("-" * 80)
        frequencies, counts = self.calculate_character_frequencies(sequences)

        total_chars = sum(counts.values())
        print(f"Total characters generated: {total_chars}")
        print(f"\n{'Character':<12} {'Count':<10} {'Frequency':<12} {'Bar'}")
        print("-" * 80)

        for char in sorted(frequencies.keys()):
            freq = frequencies[char]
            count = counts[char]
            bar = '#' * int(freq * 50)
            print(f"{char:<12} {count:<10} {freq:>10.4f}  {bar}")

        # Calculate and display empirical transitions
        print("\n" + "=" * 80)
        print("\nEmpirical Transition Probabilities (from generated sequences):")
        print("-" * 80)

        empirical_matrix = self.calculate_empirical_transitions(sequences)

        for state in sorted(empirical_matrix.keys()):
            print(f"\nFrom '{state}':")
            for next_state, prob in sorted(empirical_matrix[state].items()):
                print(f"  -> '{next_state}': {prob:.4f}")

        # Compare matrices
        print("\n" + "=" * 80)
        print("\nValidation: Empirical vs. Expected Transition Probabilities")
        print("-" * 80)

        differences = self.compare_matrices(empirical_matrix)

        if show_full_comparison:
            for state in sorted(differences.keys()):
                print(f"\nState '{state}':")
                print(f"  {'Next State':<12} {'Expected':<12} {'Observed':<12} {'Difference':<12}")
                print("  " + "-" * 60)
                for next_state in sorted(differences[state].keys()):
                    stats = differences[state][next_state]
                    print(f"  {next_state:<12} {stats['expected']:>10.4f}  {stats['observed']:>10.4f}  {stats['difference']:>10.4f}")

        # Calculate average difference
        all_diffs = []
        for state_diffs in differences.values():
            for stats in state_diffs.values():
                all_diffs.append(stats['difference'])

        if all_diffs:
            avg_diff = np.mean(all_diffs)
            max_diff = np.max(all_diffs)
            print(f"\nAverage absolute difference: {avg_diff:.4f}")
            print(f"Maximum absolute difference: {max_diff:.4f}")

            if avg_diff < 0.05:
                print("\n[OK] VALIDATION PASSED: Generated sequences closely match the transition matrix!")
            elif avg_diff < 0.1:
                print("\n[OK] VALIDATION GOOD: Generated sequences reasonably match the transition matrix.")
            else:
                print("\n[!] Note: Larger differences may occur with small sample sizes.")

        print("\n" + "=" * 80)


def load_transition_matrix(filename):
    """
    Load transition matrix from JSON file.

    Args:
        filename: Path to JSON file

    Returns:
        Transition matrix dictionary
    """
    with open(filename, 'r') as f:
        matrix = json.load(f)
    return matrix


def main():
    """Main application entry point."""
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Check if matrix file is provided
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        if not os.path.isabs(input_file):
            if os.path.exists(input_file):
                input_file = os.path.abspath(input_file)
            else:
                input_file = os.path.join(script_dir, input_file)

        print(f"Loading transition matrix from {input_file}...")
        transition_matrix = load_transition_matrix(input_file)
    else:
        print("No input file provided. Using default example matrix...")
        print("Usage: python markov_generator.py [transition_matrix.json]\n")

        # Example DNA-like transition matrix
        transition_matrix = {
            "A": {"A": 0.1, "C": 0.3, "G": 0.4, "T": 0.2},
            "C": {"A": 0.2, "C": 0.2, "G": 0.3, "T": 0.3},
            "G": {"A": 0.3, "C": 0.3, "G": 0.1, "T": 0.3},
            "T": {"A": 0.4, "C": 0.2, "G": 0.2, "T": 0.2}
        }

    # Create generator
    generator = MarkovSequenceGenerator(transition_matrix)

    # Generate sequences of varying lengths
    print("\nGenerating sequences...")

    # Define sequence lengths
    lengths = [10, 10, 20, 20, 50, 50, 50, 100, 100, 200]
    num_sequences = len(lengths)

    # Generate sequences
    sequences = generator.generate_multiple_sequences(num_sequences, lengths)

    # Display comprehensive results
    generator.display_results(sequences, show_full_comparison=True)

    # Save generated sequences to file
    output_file = os.path.join(script_dir, 'generated_sequences.txt')
    with open(output_file, 'w') as f:
        for i, seq in enumerate(sequences, 1):
            f.write(f"Sequence {i} (length={len(seq)}):\n{seq}\n\n")

    print(f"\nGenerated sequences saved to: {output_file}")


if __name__ == "__main__":
    main()
