"""
English Text Word Transition Matrix Calculator
Computes transition probabilities between words in English text.
Each word is represented by a single ASCII symbol for simplified matrix representation.
"""

import json
import os
import sys
import re
from collections import defaultdict


class WordTransitionCalculator:
    """
    Calculates transition probabilities between words in English text.
    Maps each unique word to a single ASCII symbol for matrix representation.
    """

    def __init__(self, text):
        """
        Initialize with English text.

        Args:
            text: String containing English text with spaces, punctuation, etc.
        """
        self.original_text = text
        self.words = []
        self.word_to_symbol = {}
        self.symbol_to_word = {}
        self.transition_counts = defaultdict(lambda: defaultdict(int))
        self.transition_matrix = None
        self.unique_words = []

        # Process the text
        self._extract_words()
        self._create_word_mappings()

    def _extract_words(self):
        """Extract words from text, preserving order."""
        # Split by whitespace and punctuation, but keep the words
        # Convert to lowercase for consistency
        words = re.findall(r'\b\w+\b', self.original_text.lower())
        self.words = words

        if len(self.words) < 2:
            raise ValueError("Text must contain at least 2 words")

    def _create_word_mappings(self):
        """Create bidirectional mapping between words and ASCII symbols."""
        # Get unique words in order of first appearance
        seen = set()
        for word in self.words:
            if word not in seen:
                self.unique_words.append(word)
                seen.add(word)

        # Map each word to a unique ASCII symbol
        # Start with printable ASCII characters (33-126, excluding some special chars)
        # Use extended ASCII if needed
        ascii_start = 65  # Start with 'A'

        for i, word in enumerate(self.unique_words):
            # Use letters, then numbers, then extended characters
            if i < 26:  # A-Z
                symbol = chr(65 + i)
            elif i < 52:  # a-z
                symbol = chr(97 + (i - 26))
            elif i < 62:  # 0-9
                symbol = chr(48 + (i - 52))
            elif i < 94:  # Special printable ASCII
                symbol = chr(33 + (i - 62))
            else:  # Extended ASCII
                symbol = chr(161 + (i - 94))

            self.word_to_symbol[word] = symbol
            self.symbol_to_word[symbol] = word

    def calculate_transitions(self):
        """
        Count transitions between consecutive words.
        """
        # Count all transitions
        for i in range(len(self.words) - 1):
            current_word = self.words[i]
            next_word = self.words[i + 1]

            current_symbol = self.word_to_symbol[current_word]
            next_symbol = self.word_to_symbol[next_word]

            self.transition_counts[current_symbol][next_symbol] += 1

        # Calculate probabilities
        self.transition_matrix = []
        symbols = [self.word_to_symbol[word] for word in self.unique_words]

        for from_symbol in symbols:
            row = []
            total_transitions = sum(self.transition_counts[from_symbol].values())

            for to_symbol in symbols:
                if total_transitions > 0:
                    probability = self.transition_counts[from_symbol][to_symbol] / total_transitions
                else:
                    probability = 0.0
                row.append(round(probability, 4))

            self.transition_matrix.append(row)

        return self.transition_matrix

    def get_transition_counts(self):
        """Return the raw transition counts with word labels."""
        counts_dict = {}
        for from_symbol, transitions in self.transition_counts.items():
            from_word = self.symbol_to_word[from_symbol]
            counts_dict[from_word] = {}
            for to_symbol, count in transitions.items():
                to_word = self.symbol_to_word[to_symbol]
                counts_dict[from_word][to_word] = count
        return counts_dict

    def display_results(self):
        """Display the transition matrix and statistics."""
        print("=" * 80)
        print("ENGLISH TEXT WORD TRANSITION MATRIX CALCULATOR")
        print("=" * 80)
        print(f"\nOriginal Text ({len(self.original_text)} characters):")
        print(self.original_text[:200] + ("..." if len(self.original_text) > 200 else ""))
        print(f"\nExtracted {len(self.words)} words, {len(self.unique_words)} unique")
        print("\n" + "=" * 80)

        print("\nWord to Symbol Mapping:")
        print("-" * 80)
        for word in self.unique_words:
            symbol = self.word_to_symbol[word]
            count = self.words.count(word)
            print(f"  '{symbol}' = '{word}' (appears {count} times)")

        print("\n" + "=" * 80)
        print("\nTransition Counts (word-based):")
        print("-" * 80)

        counts = self.get_transition_counts()
        for from_word in self.unique_words:
            if from_word in counts and counts[from_word]:
                print(f"\nFrom '{from_word}':")
                for to_word, count in counts[from_word].items():
                    if count > 0:  # Only show non-zero transitions
                        print(f"  -> '{to_word}': {count}")

        print("\n" + "=" * 80)
        print("\nTransition Probability Matrix (symbol-based):")
        print("-" * 80)

        symbols = [self.word_to_symbol[word] for word in self.unique_words]

        # Header
        print(f"{'From/To':<10}", end="")
        for symbol in symbols[:20]:  # Show first 20 symbols
            print(f"{symbol:>8}", end="")
        if len(symbols) > 20:
            print("  ...", end="")
        print()
        print("-" * 80)

        # Matrix rows (show first 20 rows)
        for i, from_symbol in enumerate(symbols[:20]):
            print(f"{from_symbol:<10}", end="")
            for j, to_symbol in enumerate(symbols[:20]):
                prob = self.transition_matrix[i][j]
                if prob > 0:
                    print(f"{prob:>8.4f}", end="")
                else:
                    print(f"{'':>8}", end="")
            if len(symbols) > 20:
                print("  ...", end="")
            print()

        if len(symbols) > 20:
            print("...")

        print("\n" + "=" * 80)
        print(f"\nNote: Full {len(symbols)}×{len(symbols)} matrix saved to JSON file")
        print("=" * 80)

    def save_to_json(self, output_file):
        """
        Save transition matrix and metadata to JSON file.

        Args:
            output_file: Path to output JSON file
        """
        symbols = [self.word_to_symbol[word] for word in self.unique_words]

        output_data = {
            "original_text": self.original_text,
            "text_length": len(self.original_text),
            "total_words": len(self.words),
            "unique_words_count": len(self.unique_words),
            "word_to_symbol_mapping": self.word_to_symbol,
            "symbol_to_word_mapping": self.symbol_to_word,
            "unique_words": self.unique_words,
            "symbols": symbols,
            "transition_counts": self.get_transition_counts(),
            "transition_matrix": self.transition_matrix,
            "matrix_description": f"Matrix size: {len(symbols)}×{len(symbols)}. Rows and columns ordered by unique_words list."
        }

        with open(output_file, 'w') as f:
            json.dump(output_data, f, indent=2)

        print(f"\nTransition matrix saved to {output_file}")


def load_text_from_file(filename):
    """Load English text from a text file."""
    with open(filename, 'r', encoding='utf-8') as f:
        text = f.read()
    return text


def generate_sample_text():
    """Generate sample English text for testing."""
    return """The quick brown fox jumps over the lazy dog. The dog was sleeping under a tree.
    The tree provided shade from the hot sun. The sun shone brightly in the clear sky.
    The sky was blue and beautiful. The beautiful day made everyone happy. The happy children
    played in the park. The park had many trees and flowers. The flowers bloomed in spring.
    The spring brought new life to the garden. The garden was full of colors and scents.
    The scents attracted butterflies and bees. The bees collected nectar from the flowers."""


def main():
    """Main application entry point."""
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Check if text file is provided
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        if not os.path.isabs(input_file):
            if os.path.exists(input_file):
                input_file = os.path.abspath(input_file)
            else:
                input_file = os.path.join(script_dir, input_file)

        print(f"Loading English text from {input_file}...")
        text = load_text_from_file(input_file)
    else:
        # Use sample text if no input provided
        print("No input file provided. Using sample English text...")
        print("Usage: python word_transition.py [text_file.txt]\n")
        text = generate_sample_text()

    # Create calculator and compute transitions
    calculator = WordTransitionCalculator(text)
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
        output_file = os.path.join(script_dir, 'word_transition_matrix.json')

    # Save results
    calculator.save_to_json(output_file)


if __name__ == "__main__":
    main()
