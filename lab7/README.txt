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
The project includes analysis of both a human insulin gene sequence and 10
influenza virus genomes.


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

Exercise 3: Influenza Genomes Analysis
  Download and analyze 10 influenza virus genomes from NCBI, generating
  frequency plots for each genome and comparative analysis.


FILES
----------

Main Analysis Files:
- dna_sequence.txt - Human insulin gene sequence from NCBI (1,671 bp)
- dna_repetition_detector.py - Main application for single sequence analysis
- repetition_analysis.txt - Detailed analysis results (generated)
- repetition_frequency_plots.png - Visualization of repetition frequencies

Influenza Genomes Analysis:
- analyze_influenza_genomes.py - Script for analyzing multiple genomes
- influenza_genomes/ - Directory containing 10 influenza genome sequences
- influenza_analysis_results/ - Directory with all generated plots


HOW TO RUN
----------

Single Sequence Analysis (Human Insulin Gene):
  cd lab7
  python dna_repetition_detector.py

Multiple Genomes Analysis (10 Influenza Viruses):
  cd lab7
  python analyze_influenza_genomes.py


FEATURES
----------

Core Functionality:
  - Reads DNA sequences in FASTA format
  - Detects all repetitive patterns from 6 to 10 base pairs
  - Identifies positions of each repetition in the sequence
  - Filters patterns that appear at least 2 times
  - Generates frequency plots and visualizations
  - Supports batch analysis of multiple genomes
  - Creates comparative analysis plots across genomes

Analysis Output:
  1. Console output: Summary with top 50 most frequent repetitions
  2. File output: Complete detailed analysis saved to repetition_analysis.txt
  3. Visual output: Frequency plots saved as PNG images
  4. Comparative plots: Summary comparisons across multiple genomes

Statistics Provided:
  - Total unique repetitive patterns found
  - Pattern frequencies (number of occurrences)
  - Exact positions of each repetition
  - Statistics breakdown by pattern length (6-10 bp)
  - Visual representation of frequency distribution
  - Cross-genome comparison metrics


RESULTS SUMMARY - HUMAN INSULIN GENE
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


RESULTS SUMMARY - INFLUENZA GENOMES
----------

Total Genomes Analyzed: 10 influenza virus genomes from NCBI

Genome Sources (NCBI Accession Numbers):
  1. NC_007373.1 - A/New York/392/2004(H3N2)
  2. CY121496.1 - A/California/07/2009(H1N1)
  3. KF021598.1 - A/Shanghai/02/2013(H7N9)
  4. CY033625.1 - A/Puerto Rico/8/1934(H1N1)
  5. CY125945.1 - A/Texas/50/2012(H3N2)
  6. CY147711.1 - A/Wisconsin/67/2005(H3N2)
  7. KJ942680.1 - A/Anhui/1/2013(H7N9)
  8. MF278872.1 - A/Michigan/45/2015(H1N1)
  9. MK629896.1 - A/Brisbane/02/2018(H1N1)
  10. MW626058.1 - A/Sydney/5/2021(H3N2)

Key Findings:
  - H3N2/H7N9 genomes: 1295 bp, ~383 unique patterns each
  - H1N1 genomes: 923-924 bp, ~229 unique patterns each
  - Most frequent patterns in H3N2: AACAGA (4 occurrences)
  - Most frequent patterns in H1N1: ACCAAA (6 occurrences)

Generated Outputs:
  - 10 individual frequency plots (one per genome)
  - 1 comparative summary plot showing cross-genome analysis
  - All plots include:
      * Top 20 most frequent patterns
      * Frequency distribution histogram
      * Statistics by pattern length


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

Human Insulin Gene:
- Accession: NM_001276760.2
- Organism: Homo sapiens (human)
- Gene: INS (insulin)
- Transcript: variant 1
- Type: mRNA sequence

Influenza Genomes:
- Source: NCBI GenBank
- Segment: 4 (hemagglutinin protein)
- Subtypes: H1N1, H3N2, H7N9
- Time span: 1934-2021
- All sequences are authentic NCBI reference sequences


NOTES
----------
- The algorithm is case-insensitive (converts all sequences to uppercase)
- Only standard nucleotides (A, T, G, C) are considered
- Repetitions can overlap in the sequence
- The detailed results file contains all patterns sorted by frequency
- Plots are automatically generated and saved as PNG files
- All genome sequences are downloaded from NCBI and verifiable by accession


VISUALIZATION
----------

Single Sequence Analysis Plots:
  1. Top 20 Most Frequent Patterns - Bar chart showing the most common patterns
  2. Frequency Distribution - Histogram of pattern occurrence frequencies
  3. Statistics by Pattern Length - Comparison of unique patterns vs occurrences

Multiple Genomes Analysis Plots:
  Individual genome plots (10 total):
    - Top 20 patterns for each genome
    - Frequency distribution for each genome
    - Pattern length statistics for each genome

  Comparative summary plot:
    - Total unique patterns per genome
    - Genome sequence lengths comparison
    - Pattern distribution by length across all genomes
    - Highest pattern repetition counts comparison


PROJECT STRUCTURE
----------
lab7/
  |-- dna_sequence.txt
  |-- dna_repetition_detector.py
  |-- repetition_analysis.txt
  |-- repetition_frequency_plots.png
  |-- analyze_influenza_genomes.py
  |-- README.txt
  |-- influenza_genomes/
  |     |-- genome_1.txt
  |     |-- genome_2.txt
  |     |-- ... (genome_3 through genome_10)
  |-- influenza_analysis_results/
        |-- Genome_1_frequency_plot.png
        |-- Genome_2_frequency_plot.png
        |-- ... (plots for all 10 genomes)
        |-- genomes_comparison_summary.png

================================================================================
