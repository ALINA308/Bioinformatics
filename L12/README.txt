================================================================================
BIOINFORMATICS LABORATORY 12 - SOLUTIONS
================================================================================

Student: Bilciurescu Elena-Alina
Group: 1241EA


================================================================================
EXERCISE 1: Matrix Prediction (Discrete Step Predictions)
================================================================================

Location: ex1/

Problem:
Implement a software application that performs predictions using a square
matrix of arbitrary size and an initial vector. The application makes
predictions over 5 discrete steps.

Files:
  - predictor.py          : Main application for matrix-based predictions
  - example_input.json    : Sample 3x3 matrix with initial vector
  - example_markov.json   : Markov chain example
  - test_2x2.json         : 2x2 test case
  - test_4x4.json         : 4x4 test case
  - results.json          : Output file (auto-generated)
  - README.md             : Complete documentation

Implementation:
  - Accepts arbitrary-sized square matrices
  - Performs exactly 5 discrete prediction steps
  - Uses matrix multiplication: State(t+1) = Matrix × State(t)
  - Validates input dimensions
  - Outputs results to JSON format

Usage:
  cd ex1
  python predictor.py example_input.json

Status: COMPLETED AND TESTED


================================================================================
EXERCISE 2: DNA Transition Matrix Calculator
================================================================================

Location: ex2/

Problem:
Use a random DNA sequence of approximately 50 letters. Calculate the transition
probabilities between nucleotides (A, C, G, T). Output the transition matrix
stored as a JSON file.

Files:
  - dna_transition.py       : DNA sequence analyzer
  - sample_dna.txt          : 49-nucleotide sample sequence
  - transition_matrix.json  : Output 4x4 transition matrix (auto-generated)
  - README.md               : Complete documentation

Implementation:
  - Processes DNA sequences (A, C, G, T nucleotides)
  - Counts transitions between consecutive nucleotides
  - Calculates probability matrix (4×4)
  - Validates DNA sequence format
  - Can generate random sequences for testing

Usage:
  cd ex2
  python dna_transition.py sample_dna.txt

Output Format:
  - Transition matrix: P(nucleotide_i -> nucleotide_j)
  - Row-stochastic matrix (rows sum to 1.0)
  - JSON format with detailed statistics

Status: COMPLETED AND TESTED


================================================================================
EXERCISE 3: English Text Word Transition Matrix
================================================================================

Location: ex3/

Problem:
Use a random English text of approximately 300 letters (including spaces,
punctuation, and everything). Compute the transition probabilities between
consecutive words. Store the transition matrix as a JSON file. Each unique
word is represented by a single ASCII symbol.

Files:
  - word_transition.py          : Word transition analyzer
  - sample_text.txt             : 420-character English text sample
  - word_transition_matrix.json : Output matrix with word mappings (auto-generated)
  - README.md                   : Complete documentation

Implementation:
  - Extracts words using regex pattern matching
  - Maps each unique word to ASCII symbol (A-Z, a-z, 0-9, ...)
  - Calculates N×N transition matrix (N = unique words)
  - Case-insensitive processing
  - Complete word-to-symbol bidirectional mapping

Usage:
  cd ex3
  python word_transition.py sample_text.txt

Output Format:
  - Word-to-symbol mapping
  - Symbol-to-word mapping
  - Transition matrix with probabilities
  - Transition counts for validation

Status: COMPLETED AND TESTED


================================================================================
EXERCISE 4: Markov Chain Sequence Generator
================================================================================

Location: ex4/

Problem:
Implement a Markov chain-based sequence generator that reads a transition
matrix from a JSON file and generates new sequences based on learned transition
probabilities. The application must validate the matrix, generate multiple
sequences of varying lengths, display statistics, and compare empirical vs.
expected transition probabilities.

Files:
  - markov_generator.py       : Markov chain sequence generator
  - dna_matrix.json           : DNA 4-state transition matrix (A, C, G, T)
  - protein_matrix.json       : Protein 4-state matrix (A, B, P, E)
  - generated_sequences.txt   : Generated sequences output (auto-generated)
  - README.md                 : Complete documentation

Implementation:
  - Loads and validates transition matrices from JSON
  - Validates probabilities (sum to 1.0, non-negative)
  - Generates sequences using np.random.choice() with probabilities
  - Creates 10 sequences of varying lengths (10, 20, 50, 100, 200)
  - Calculates character frequency distributions
  - Computes empirical transition probabilities
  - Validates by comparing empirical vs. expected probabilities

Usage:
  cd ex4
  python markov_generator.py dna_matrix.json

Features:
  - Random or specified starting states
  - Statistical validation with average/max differences
  - Visual frequency bars
  - Complete comparison tables

Status: COMPLETED AND TESTED


================================================================================
TECHNICAL DETAILS
================================================================================

Programming Language: Python 3.11
External Dependencies:
  - NumPy (for matrix operations and probability sampling)

All applications:
  - Include comprehensive error handling
  - Validate input data
  - Provide detailed output statistics
  - Save results to JSON or text files
  - Include complete documentation in README.md files

Testing:
  - All exercises tested with sample data
  - Validation performed on outputs
  - Edge cases handled appropriately


================================================================================
RUNNING THE APPLICATIONS
================================================================================

Note for Windows PowerShell Users:
PowerShell does not support '&&' as a command separator. Use semicolon ';'
instead, or run commands separately.

Example (PowerShell):
  cd ex1 ; python predictor.py example_input.json

Or run separately:
  cd ex1
  python predictor.py example_input.json



These implementations are applicable to:

1. Markov Chain Modeling: State transition predictions
2. DNA Sequence Analysis: Nucleotide transition patterns
3. Natural Language Processing: Word co-occurrence patterns
4. Sequence Generation: Synthetic data creation
5. Gene Prediction: Identifying coding regions
6. Protein Modeling: Amino acid transitions
7. Evolutionary Analysis: Sequence evolution simulation


