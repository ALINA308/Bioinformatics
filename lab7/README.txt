================================================================================
Lab 7: DNA Sequence Repetition Detection
================================================================================

AUTHOR:
----------
Bilciurescu Elena Alina 1241EA


OVERVIEW
----------
This project implements a software application that detects repetitive patterns
in DNA sequences from NCBI (National Center for Biotechnology Information).


EXERCISES IMPLEMENTED
----------

Exercise 1: DNA Sequence from NCBI
  - Source: NCBI Reference Sequence NM_001276760.2
  - Gene: Homo sapiens insulin (INS), transcript variant 1, mRNA
  - Length: 1,671 base pairs (nucleotides)
  - File: dna_sequence.txt

Exercise 2: Repetition Detection Software
  The application detects all repetitions between 6 and 10 base pairs in the
  DNA sequence.


FILES
----------
- dna_sequence.txt - The DNA sequence from NCBI (1,671 bp)
- dna_repetition_detector.py - Main application for detecting repetitions
- repetition_analysis.txt - Detailed analysis results (generated)
- repetition_frequency_plots.png - Visualization of repetition frequencies


HOW TO RUN
----------
  cd lab7
  python dna_repetition_detector.py


FEATURES
----------

Core Functionality:
  - Reads DNA sequences in FASTA format
  - Detects all repetitive patterns from 6 to 10 base pairs
  - Identifies positions of each repetition in the sequence
  - Filters patterns that appear at least 2 times
  - Generates frequency plots and visualizations

Analysis Output:
  1. Console output: Summary with top 50 most frequent repetitions
  2. File output: Complete detailed analysis saved to repetition_analysis.txt
  3. Visual output: Frequency plots saved as PNG image

Statistics Provided:
  - Total unique repetitive patterns found
  - Pattern frequencies (number of occurrences)
  - Exact positions of each repetition
  - Statistics breakdown by pattern length (6-10 bp)
  - Visual representation of frequency distribution


RESULTS SUMMARY
----------

Key Findings:
  - Total unique patterns: 1,299
  - Most frequent pattern: GCGGGG (23 occurrences)
  - Pattern length distribution:
      * 6 bp: 348 unique patterns (1,223 total occurrences)
      * 7 bp: 326 unique patterns (924 total occurrences)
      * 8 bp: 257 unique patterns (648 total occurrences)
      * 9 bp: 200 unique patterns (467 total occurrences)
      * 10 bp: 168 unique patterns (371 total occurrences)

Top 5 Most Frequent Repetitions:
  1. GCGGGG - 23 occurrences
  2. GGCGGG - 20 occurrences
  3. CCCGCC - 17 occurrences
  4. GGGCGG - 17 occurrences
  5. CGGGGG - 16 occurrences


ALGORITHM DETAILS
----------

Repetition Detection Algorithm:

  1. Input Processing: Reads DNA sequence from FASTA file, removing headers

  2. Pattern Scanning: For each pattern length (6-10 bp):
     - Slides a window through the entire sequence
     - Extracts each possible substring of that length
     - Searches for all occurrences of that substring

  3. Filtering: Keeps only patterns that appear 2+ times

  4. Analysis: Sorts patterns by frequency and generates statistics

  5. Output: Displays results and saves detailed report

Time Complexity:
  - Overall: O(n × m × L) where:
      * n = sequence length
      * m = sum of pattern lengths (6+7+8+9+10 = 40)
      * L = average pattern length
  - Space Complexity: O(k) where k = number of unique patterns found


REQUIREMENTS
----------
- Python 3.x
- matplotlib (for plotting)


DNA SEQUENCE INFORMATION
----------
- Accession: NM_001276760.2
- Organism: Homo sapiens (human)
- Gene: INS (insulin)
- Transcript: variant 1
- Type: mRNA sequence


NOTES
----------
- The algorithm is case-insensitive (converts all sequences to uppercase)
- Only standard nucleotides (A, T, G, C) are considered
- Repetitions can overlap in the sequence
- The detailed results file contains all patterns sorted by frequency
- Plots are automatically generated and saved as PNG files


VISUALIZATION
----------
The application generates three plots:
  1. Top 20 Most Frequent Patterns - Bar chart showing the most common patterns
  2. Frequency Distribution - Histogram of pattern occurrence frequencies
  3. Statistics by Pattern Length - Comparison of unique patterns vs occurrences

================================================================================
