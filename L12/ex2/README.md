# DNA Sequence Transition Matrix Calculator

This application computes transition probabilities between nucleotides in a DNA sequence. It analyzes consecutive nucleotide pairs to build a 4×4 transition matrix showing the probability of transitioning from one nucleotide (A, C, G, T) to another.

## Problem Description

Given a DNA sequence of approximately 50 letters, calculate the transition probabilities between consecutive nucleotides and output the transition matrix as a JSON file.

## Features

- Processes DNA sequences containing A, C, G, T nucleotides
- Counts all transitions between consecutive nucleotides
- Calculates transition probabilities
- Outputs a 4×4 transition matrix (A→C→G→T order)
- Saves results to JSON file with detailed statistics
- Validates input sequences
- Can generate random DNA sequences for testing

## Requirements

- Python 3.6+
- No external dependencies (uses only Python standard library)

## Usage

**IMPORTANT:** Run from within the `ex2` folder:

```bash
cd ex2
python dna_transition.py [sequence_file.txt] [output_file.json]
```

### Basic Usage (Random Sequence)

Generate and analyze a random 50-nucleotide DNA sequence:
```bash
python dna_transition.py
```

This creates `transition_matrix.json` automatically.

### Using Input File

Analyze a DNA sequence from a text file:
```bash
python dna_transition.py sample_dna.txt
```

### Custom Output File

Specify a custom output filename:
```bash
python dna_transition.py sample_dna.txt my_results.json
```

## Input Format

Create a text file containing your DNA sequence:

```
ACGTACGTTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
```

- Valid nucleotides: A, C, G, T (case-insensitive)
- Minimum length: 2 nucleotides
- Recommended: ~50 nucleotides for meaningful statistics

## Output Format

The `transition_matrix.json` file contains:

```json
{
  "dna_sequence": "ACGT...",
  "sequence_length": 50,
  "nucleotides": ["A", "C", "G", "T"],
  "transition_counts": {
    "A": {"A": 0, "C": 5, "G": 3, "T": 2},
    "C": {"A": 1, "C": 0, "G": 8, "T": 1},
    "G": {"A": 2, "C": 6, "G": 0, "T": 4},
    "T": {"A": 7, "C": 2, "G": 3, "T": 0}
  },
  "transition_matrix": [
    [0.0, 0.5, 0.3, 0.2],
    [0.1, 0.0, 0.8, 0.1],
    [0.17, 0.5, 0.0, 0.33],
    [0.58, 0.17, 0.25, 0.0]
  ],
  "matrix_description": "Rows and columns ordered as: A, C, G, T"
}
```

### Matrix Interpretation

- **Rows**: Source nucleotide (from)
- **Columns**: Target nucleotide (to)
- **Order**: A, C, G, T
- **Example**: `transition_matrix[0][1]` = probability of A→C transition

## Example Files

### sample_dna.txt
A 50-nucleotide DNA sequence for testing.

## Examples

All commands assume you are in the `ex2` folder:

### Analyze sample DNA:
```bash
python dna_transition.py sample_dna.txt
```

### Generate random sequence:
```bash
python dna_transition.py
```

## Output Display

The program displays:
1. Original DNA sequence
2. Transition counts table
3. Transition probability matrix
4. Confirmation of JSON file creation

## Bioinformatics Applications

DNA transition matrices are used in:

1. **Sequence Analysis**: Understanding nucleotide composition patterns
2. **Markov Models**: Building probabilistic models of DNA sequences
3. **Gene Finding**: Identifying coding vs. non-coding regions
4. **Mutation Analysis**: Studying nucleotide substitution patterns
5. **Sequence Generation**: Creating synthetic DNA sequences

## Implementation Details

- Uses dictionary-based counting for transitions
- Probabilities calculated as: P(from→to) = count(from→to) / total_transitions(from)
- Handles edge cases (no transitions from a nucleotide)
- Rounds probabilities to 4 decimal places
- All results saved to ex2 folder automatically

## Quick Start

1. Navigate to ex2 folder:
   ```bash
   cd C:\Users\alina\Desktop\L12\ex2
   ```

2. Run with sample DNA:
   ```bash
   python dna_transition.py sample_dna.txt
   ```

3. Check `transition_matrix.json` for results

## Mathematical Model

For each nucleotide pair (n_i, n_{i+1}) in the sequence:
- Count transitions: count[n_i][n_{i+1}]++
- Calculate probability: P(n_i → n_{i+1}) = count[n_i][n_{i+1}] / Σ_j count[n_i][n_j]

The resulting 4×4 matrix is row-stochastic (each row sums to 1.0 or 0.0 if no transitions).
