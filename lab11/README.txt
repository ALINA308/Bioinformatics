================================================================================
                         LABORATORY 11 - MARKOV MODELS
================================================================================

Student: Bilciurescu Elena-Alina
Group: 1241EA


================================================================================
EXERCISE 1: CpG ISLAND DETECTION USING MARKOV MODELS
================================================================================

PROBLEM DESCRIPTION:
-------------------
Identify if a DNA sequence is a CpG island based on transition probabilities
between nucleotides using a Position Weight Matrix (PWM) approach.

TRAINING DATA:
- Positive Sequence (Island): ATCGATTCGATATCATACACGTAT
- Negative Sequence (Non-Island): CTCGACTAGTATGAAGTCCACGCTTG
- Test Sequence: CAGGTTGGAAACGTAA

SOLUTION APPROACH:
-----------------
1. TRANSITION COUNTING:
   - Created two 4x4 transition frequency matrices (A, C, G, T)
   - Counted every occurrence where nucleotide i is followed by nucleotide j
   - Applied pseudocount (epsilon=1e-5) to avoid division by zero

2. PROBABILITY CALCULATION:
   - Converted counts to probabilities by normalizing each row
   - Each cell [i,j] represents P(nucleotide_j | nucleotide_i)

3. LOG-LIKELIHOOD MATRIX:
   - Formula: log2(Prob_positive / Prob_negative)
   - Positive values indicate CpG island transitions
   - Negative values indicate non-island transitions

4. SEQUENCE SCORING:
   - Calculated total score by summing log-likelihood values for each transition
   - Decision rule: Score > 0 → CpG Island, Score < 0 → Non-Island

RESULTS:
--------
The algorithm successfully classifies DNA sequences based on their transition
patterns, distinguishing between CpG island and non-island regions.

FILE: ex1.py


================================================================================
EXERCISE 2: PLAGIARISM DETECTION USING STYLOMETRIC ANALYSIS
================================================================================

PROBLEM DESCRIPTION:
-------------------
Courtroom scenario: Detect plagiarism in a suspect text by analyzing writing
style using Markov models trained on two Romanian poets:
- Mihai Eminescu (classic, romantic style)
- Nichita Stanescu (modern, abstract style)

SOLUTION APPROACH:
-----------------
1. TRAINING DATA PREPARATION:
   - Collected poetry samples from both authors
   - Eminescu: 151 words with nature imagery (luna, noapte, florile)
   - Stanescu: 177 words with existential themes (singur, gol, nu stiu)

2. TEXT PREPROCESSING:
   - Converted text to lowercase
   - Removed special characters and punctuation
   - Split into word tokens

3. TRANSITION MATRIX CONSTRUCTION:
   - Built word-to-word transition probability matrices
   - Formula: P(word_j | word_i) for each author
   - Applied epsilon=1e-5 pseudocount for unseen transitions

4. LOG-LIKELIHOOD MATRIX:
   - Combined the two probability matrices
   - Formula: log2(P_eminescu / P_stanescu)
   - Positive scores → Eminescu-like style
   - Negative scores → Stanescu-like style

5. SLIDING WINDOW ANALYSIS:
   - Window size: 10 words
   - Step size: 3 words
   - Calculated average log-likelihood score per window

6. CLASSIFICATION:
   - Score > 0.5  → EMINESCU style
   - Score < -0.5 → STANESCU style
   - -0.5 ≤ Score ≤ 0.5 → UNCERTAIN/MIXED

RESULTS:
--------
Common vocabulary: 7 words (si, in, este, sunt, totul, nu, foarte)

VERDICT SUMMARY:
- Windows attributed to Eminescu: 6
- Windows attributed to Stanescu: 3
- Uncertain/Mixed windows: 20

CONCLUSION: The accused text contains passages from BOTH authors.
Evidence of PLAGIARISM detected!

The algorithm successfully identified:
- Eminescu-style passages (romantic/nature imagery)
- Stanescu-style passages (existential/abstract themes)
- Mixed style passages

FILE: ex2.py


================================================================================
TECHNICAL IMPLEMENTATION DETAILS
================================================================================

LIBRARIES USED:
- numpy: Numerical computations and logarithms
- collections.defaultdict: Efficient transition counting
- re: Text preprocessing and cleaning

KEY ALGORITHMS:
1. Position Weight Matrix (PWM) for sequence analysis
2. First-order Markov Model for transition probabilities
3. Log-likelihood ratio test for classification
4. Sliding window technique for stylometric analysis

MATHEMATICAL FOUNDATIONS:
- Probability theory and conditional probabilities
- Logarithmic transformations to avoid numerical underflow
- Statistical classification using threshold-based decision rules


