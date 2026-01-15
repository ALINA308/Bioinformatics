# English Text Word Transition Matrix Calculator

This application computes transition probabilities between consecutive words in English text. Each unique word is mapped to a single ASCII symbol for simplified matrix representation, making it easier to visualize and work with the transition matrix.

## Problem Description

Given an English text of approximately 300 letters (including spaces, punctuation, and everything), compute the transition probabilities between consecutive words. Store the transition matrix as a JSON file. Each unique word is represented by a single ASCII symbol.

## Features

- Processes English text with spaces, punctuation, and all characters
- Extracts words using regex pattern matching
- Maps each unique word to a unique ASCII symbol (A-Z, a-z, 0-9, then special characters)
- Counts transitions between consecutive words
- Calculates transition probabilities
- Outputs an N×N transition matrix (where N = number of unique words)
- Saves results to JSON file with complete word mappings
- Handles arbitrary text lengths
- Case-insensitive word processing

## Requirements

- Python 3.6+
- No external dependencies (uses only Python standard library)

## Usage

**IMPORTANT:** Run from within the `ex3` folder:

```bash
cd ex3
python word_transition.py [text_file.txt] [output_file.json]
```

### Basic Usage (Sample Text)

Analyze sample English text:
```bash
python word_transition.py
```

This creates `word_transition_matrix.json` automatically.

### Using Input File

Analyze text from a file:
```bash
python word_transition.py sample_text.txt
```

### Custom Output File

Specify a custom output filename:
```bash
python word_transition.py sample_text.txt my_results.json
```

## Input Format

Create a text file containing English text (approximately 300 characters):

```
The quick brown fox jumps over the lazy dog. The dog barked loudly at the mailman.
The mailman delivered packages to every house on the street. The street was lined with
beautiful oak trees...
```

- Include spaces, punctuation, capitalization - everything is processed
- Words are extracted using regex pattern `\b\w+\b`
- Converted to lowercase for consistency
- Minimum: 2 words required

## Output Format

The `word_transition_matrix.json` file contains:

```json
{
  "original_text": "The quick brown fox...",
  "text_length": 350,
  "total_words": 65,
  "unique_words_count": 42,
  "word_to_symbol_mapping": {
    "the": "A",
    "quick": "B",
    "brown": "C",
    "fox": "D",
    ...
  },
  "symbol_to_word_mapping": {
    "A": "the",
    "B": "quick",
    "C": "brown",
    "D": "fox",
    ...
  },
  "unique_words": ["the", "quick", "brown", "fox", ...],
  "symbols": ["A", "B", "C", "D", ...],
  "transition_counts": {
    "the": {"quick": 1, "dog": 2, "mailman": 1, ...},
    "quick": {"brown": 1},
    ...
  },
  "transition_matrix": [
    [0.0, 0.1, 0.05, 0.15, ...],  // transitions from "the"
    [0.0, 0.0, 1.0, 0.0, ...],     // transitions from "quick"
    ...
  ],
  "matrix_description": "Matrix size: 42×42. Rows and columns ordered by unique_words list."
}
```

### Matrix Interpretation

- **Rows**: Source word (from)
- **Columns**: Target word (to)
- **Order**: Same as `unique_words` list
- **Symbols**: Each word mapped to a single ASCII character
- **Example**: If "the" is mapped to 'A' and appears at index 0, `transition_matrix[0][3]` = probability of transitioning from "the" to the word at index 3

## Symbol Mapping

Words are mapped to ASCII symbols in this order:
1. **A-Z** (26 words): Uppercase letters
2. **a-z** (26 words): Lowercase letters
3. **0-9** (10 words): Digits
4. **Special ASCII** (32 words): Printable special characters
5. **Extended ASCII**: For texts with >94 unique words

## Example Files

### sample_text.txt
A ~300-character English text sample with various words and punctuation.

## Examples

All commands assume you are in the `ex3` folder:

### Analyze sample text:
```bash
python word_transition.py sample_text.txt
```

### Use default sample:
```bash
python word_transition.py
```

## Output Display

The program displays:
1. Original text (first 200 characters)
2. Word statistics (total words, unique words)
3. Word-to-symbol mapping table
4. Transition counts (word pairs)
5. Transition probability matrix (symbol-based, first 20×20 shown)
6. Confirmation of JSON file creation

## Applications

Word transition matrices are used in:

1. **Natural Language Processing**: Understanding word co-occurrence patterns
2. **Markov Text Generation**: Creating synthetic text that mimics the style
3. **Autocomplete Systems**: Predicting next word in a sequence
4. **Language Modeling**: Building probabilistic models of language
5. **Text Analysis**: Studying writing patterns and styles
6. **Predictive Keyboards**: Suggesting next words

## Implementation Details

- Uses `re.findall(r'\b\w+\b', text)` to extract words
- Converts all words to lowercase for consistency
- Dictionary-based counting for transitions
- Probabilities: P(word1→word2) = count(word1→word2) / total_transitions(word1)
- Rounds probabilities to 4 decimal places
- Matrix is row-stochastic (each row sums to 1.0 or 0.0)
- All results saved to ex3 folder automatically

## Quick Start

1. Navigate to ex3 folder:
   ```bash
   cd C:\Users\alina\Desktop\L12\ex3
   ```

2. Run with sample text:
   ```bash
   python word_transition.py sample_text.txt
   ```

3. Check `word_transition_matrix.json` for results

## Mathematical Model

For each consecutive word pair (w_i, w_{i+1}) in the text:
- Count transitions: count[w_i][w_{i+1}]++
- Calculate probability: P(w_i → w_{i+1}) = count[w_i][w_{i+1}] / Σ_j count[w_i][w_j]

The resulting N×N matrix (where N = unique words) is row-stochastic.

## Example Output

For the text "The cat and the dog. The dog and the cat."

- **Unique words**: the, cat, and, dog (N=4)
- **Symbol mapping**: A=the, B=cat, C=and, D=dog
- **Matrix**: 4×4 showing probabilities like P(the→cat)=0.5, P(the→dog)=0.5, P(cat→and)=1.0, etc.
