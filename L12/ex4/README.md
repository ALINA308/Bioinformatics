# Markov Chain-Based Sequence Generator for Bioinformatics

A Python application that generates new sequences based on learned transition probabilities using Markov chains. This tool is essential for bioinformatics applications including DNA sequence generation, protein modeling, and pattern analysis.

## Problem Description

Given a transition matrix that defines probabilities for state transitions between characters (like 'A', 'C', 'G', 'T' for DNA, or 'A', 'B', 'P', 'E' for proteins), generate new sequences that follow these learned transition patterns.

## Features

### Core Functionality
- **Load transition matrices** from JSON files
- **Validate matrices** (probabilities sum to 1.0, non-negative values)
- **Generate sequences** of any specified length
- **Multiple sequence generation** with varying lengths
- **Statistical analysis** of generated sequences
- **Validation** by comparing empirical vs. expected transition probabilities

### Advanced Features
- Random or specified starting states
- Efficient random sampling using NumPy
- Character frequency distribution analysis
- Empirical transition probability calculation
- Comprehensive comparison metrics
- Edge case handling (unknown states, invalid probabilities)

## Requirements

- Python 3.6+
- NumPy (for probability-based random sampling)

Install dependencies:
```bash
pip install numpy
```

## Usage

**IMPORTANT:** Run from within the `ex4` folder:

```bash
cd ex4
python markov_generator.py [transition_matrix.json]
```

### Basic Usage (Default Matrix)

Run with built-in DNA example:
```bash
python markov_generator.py
```

### Using Custom Transition Matrix

Analyze with DNA transition matrix:
```bash
python markov_generator.py dna_matrix.json
```

Analyze with protein transition matrix:
```bash
python markov_generator.py protein_matrix.json
```

## Input Format

Create a JSON file with transition probabilities:

```json
{
  "A": {
    "A": 0.1,
    "C": 0.3,
    "G": 0.4,
    "T": 0.2
  },
  "C": {
    "A": 0.2,
    "C": 0.2,
    "G": 0.3,
    "T": 0.3
  },
  "G": {
    "A": 0.3,
    "C": 0.3,
    "G": 0.1,
    "T": 0.3
  },
  "T": {
    "A": 0.4,
    "C": 0.2,
    "G": 0.2,
    "T": 0.2
  }
}
```

**Requirements:**
- Each state must have transition probabilities to at least one other state
- Probabilities for each state must sum to 1.0 (±0.000001 tolerance)
- All probabilities must be non-negative

## Example Files

### dna_matrix.json
4-state DNA transition matrix (A, C, G, T) with realistic transition probabilities.

### protein_matrix.json
4-state protein/abstract model (A, B, P, E) showing complex transition patterns.

## Output

### Console Output

The program displays:

1. **Original Transition Matrix**: Loaded probabilities
2. **Generated Sequences**: 10 sequences of varying lengths (10, 20, 50, 100, 200 characters)
3. **Character Frequency Statistics**: Distribution with counts and visual bars
4. **Empirical Transition Probabilities**: Calculated from generated sequences
5. **Validation Comparison**: Expected vs. observed with difference metrics

### File Output

**generated_sequences.txt**: All generated sequences saved for further analysis

## Technical Implementation

### Algorithm

1. **Matrix Validation**: Ensure probabilities are valid
2. **State Sampling**:
   - Start from random state (or specified)
   - For each step, sample next state using `np.random.choice()` with transition probabilities
3. **Sequence Generation**: Repeat for desired length
4. **Statistical Analysis**: Calculate empirical frequencies and transitions
5. **Validation**: Compare empirical vs. expected using absolute differences

### Key Methods

```python
class MarkovSequenceGenerator:
    def generate_sequence(length, start_state=None)
    def generate_multiple_sequences(num_sequences, lengths)
    def calculate_character_frequencies(sequences)
    def calculate_empirical_transitions(sequences)
    def compare_matrices(empirical_matrix)
```

## Examples

All commands assume you are in the `ex4` folder:

### Generate DNA sequences:
```bash
python markov_generator.py dna_matrix.json
```

### Generate protein sequences:
```bash
python markov_generator.py protein_matrix.json
```

### Use default example:
```bash
python markov_generator.py
```

## Validation Metrics

The program calculates:

- **Average Absolute Difference**: Mean difference between expected and observed probabilities
- **Maximum Absolute Difference**: Largest single difference
- **Validation Status**:
  - < 0.05: PASSED (excellent match)
  - < 0.10: GOOD (reasonable match)
  - ≥ 0.10: Note about sample size

## Bioinformatics Applications

This Markov chain generator is used for:

1. **DNA Sequence Simulation**: Generate synthetic DNA sequences
2. **Protein Modeling**: Model amino acid transition patterns
3. **Gene Prediction**: Learn and apply gene structure patterns
4. **Evolutionary Analysis**: Simulate sequence evolution
5. **Sequence Motif Generation**: Create sequences with specific patterns
6. **Training Data Augmentation**: Generate additional sequences for ML models
7. **Statistical Testing**: Create null model sequences for hypothesis testing

## Implementation Details

- **Efficient Sampling**: Pre-computed state and probability arrays
- **Numerical Stability**: Proper handling of floating-point probabilities
- **Validation**: Comprehensive matrix validation before generation
- **Error Handling**: Meaningful errors for invalid inputs
- **Modularity**: Well-structured OOP design
- **Documentation**: Detailed docstrings for all methods

## Example Output

For a DNA matrix with 4 states:

```
Generated Sequences:
 1. (length= 10) ACGTTAGCTA
 2. (length= 10) GTACGATCGA
 3. (length= 20) ACGATCGATCGATCGATCGA
 ...
10. (length=200) ACGTACGT...

Character Frequency Statistics:
Character    Count      Frequency    Bar
----------------------------------------------------------------
A            150        0.2500       ############
C            180        0.3000       ###############
G            120        0.2000       ##########
T            150        0.2500       ############

Validation: Empirical vs. Expected
Average absolute difference: 0.0234
Maximum absolute difference: 0.0567
[OK] VALIDATION PASSED: Generated sequences closely match the transition matrix!
```

## Quick Start

1. Navigate to ex4 folder:
   ```bash
   cd C:\Users\alina\Desktop\L12\ex4
   ```

2. Run with DNA matrix:
   ```bash
   python markov_generator.py dna_matrix.json
   ```

3. Check `generated_sequences.txt` for output

## Troubleshooting

**Error: "Probabilities sum to X, not 1.0"**
- Check that all transition probabilities for each state sum to exactly 1.0
- Small rounding errors (±0.000001) are acceptable

**Error: "Unknown starting state"**
- Ensure the specified starting state exists in the transition matrix

**Poor validation results**
- Generate more/longer sequences for better statistical accuracy
- Check that the transition matrix is correct

## Mathematical Model

For a Markov chain with states S = {s₁, s₂, ..., sₙ}:

- **Transition Matrix**: P[i,j] = P(next state = sⱼ | current state = sᵢ)
- **Generation**: At each step, sample next state according to P[current, ·]
- **Validation**: Compare empirical P'[i,j] with expected P[i,j]

The Markov property: P(Xₜ₊₁ | Xₜ, Xₜ₋₁, ..., X₁) = P(Xₜ₊₁ | Xₜ)

## PowerShell Usage Note

For Windows PowerShell users, use semicolon instead of &&:

```powershell
cd ex4 ; python markov_generator.py dna_matrix.json
```
