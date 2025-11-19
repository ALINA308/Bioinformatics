==============================================================================
TRANSPOSABLE ELEMENT DETECTION - LAB 8
==============================================================================
Bilciurescu Elena-Alina 1241 EA
OVERVIEW
--------
This project implements a bioinformatics solution for detecting transposable
elements (TEs) in DNA sequences. It consists of two main components:

1. Sequence Generator - Creates artificial DNA sequences with embedded TEs
2. TE Detector - Identifies transposable elements based on structural features


BACKGROUND
----------
Transposable elements are DNA sequences that can move within a genome. They
have characteristic structural features:

  - Terminal Inverted Repeats (TIRs): Inverted sequences at both ends
  - Target Site Duplications (TSDs): Short direct repeats flanking insertion
  - Structure: TSD - TIR - Internal Sequence - TIR_RC - TSD


FILES
-----
  generate_sequence.py           - Generates artificial DNA with 3-4 TEs
  detect_transposons.py          - Detects TEs using pattern matching
  artificial_dna.fasta           - Generated DNA sequence (output)
  transposons_ground_truth.json  - True positions of TEs (output)
  detection_results.json         - Detection results (output)


USAGE
-----

1. Generate Artificial DNA Sequence:

   python generate_sequence.py

   This creates:
   - A 200-400bp DNA sequence with 3-4 transposable elements
   - Each TE has Terminal Inverted Repeats (8-15bp)
   - Each TE has Target Site Duplications (5-12bp)


2. Detect Transposable Elements:

   python detect_transposons.py

   This analyzes the sequence and:
   - Finds inverted repeats (potential TIRs)
   - Finds direct repeats (potential TSDs)
   - Combines evidence to identify TE boundaries
   - Compares results with ground truth
   - Reports accuracy metrics


DETECTION ALGORITHM
-------------------

Step 1: Find Inverted Repeats (TIRs)
  - Scans sequence for palindromic patterns
  - Looks for sequences and their reverse complements
  - Parameters:
    * Minimum length: 8bp
    * Maximum length: 15bp
    * Maximum distance: 100bp

Step 2: Find Direct Repeats (TSDs)
  - Searches for exact sequence matches
  - Parameters:
    * Minimum length: 5bp
    * Maximum length: 12bp
    * Maximum distance: 100bp

Step 3: Combine Evidence
  - Matches TSDs that flank inverted repeat pairs
  - Determines TE boundaries
  - Assigns confidence scores:
    * High: Both TIRs and TSDs found
    * Medium: Only TIRs found


EXAMPLE OUTPUT
--------------

======================================================================
DETECTION RESULTS
----------------------------------------------------------------------
Detected 3 transposable element(s)

ID    Start    End      Length   TIR Len    TSD          Confidence
----------------------------------------------------------------------
1     10       66       57       10         CTCTGACT     high
      TIR: GGCCGAATAG ... CTATTCGGCC
      Structure: CTCTGACT|GGCCGAATAG...CTATTCGGCC|CTCTGACT

2     77       131      55       9          TGCGACAG     high
      TIR: TGACGCTTT ... AAAGCGTCA
      Structure: TGCGACAG|TGACGCTTT...AAAGCGTCA|TGCGACAG

3     142      195      54       9          TAGCAGCC     high
      TIR: GCAGTAAGG ... CCTTACTGC
      Structure: TAGCAGCC|GCAGTAAGG...CCTTACTGC|TAGCAGCC


ACCURACY METRICS
----------------
The detector achieves 100% precision and recall on generated test sequences:
  - True Positives: All inserted TEs are correctly identified
  - Precision: No false positives
  - Recall: All TEs are detected
  - Position accuracy tolerance: ±10bp


TESTING WITH OVERLAPPING TRANSPOSONS
-------------------------------------
To test edge cases (e.g., overlapping or nested transposons), modify the
parameters in generate_sequence.py:

  # Reduce spacing to create overlapping elements
  spacer_lengths = [random.randint(0, 5) for _ in range(num_transposons + 1)]


REFERENCES
----------
1. Transposable element detection from whole genome sequence data
   https://pmc.ncbi.nlm.nih.gov/articles/PMC4696183/

2. STEAK: A specific tool for transposable elements and retrovirus detection
   https://pmc.ncbi.nlm.nih.gov/articles/PMC5597868/


KEY CONCEPTS DEMONSTRATED
-------------------------
  - DNA sequence manipulation
  - Pattern matching algorithms
  - Reverse complement computation
  - Structural motif detection
  - Bioinformatics algorithm validation
  - Ground truth comparison


DEPENDENCIES
------------
  - Python 3.x
  - Standard library only (no external packages required)


==============================================================================
PART 2: BACTERIAL GENOME TRANSPOSON ANALYSIS
==============================================================================

OVERVIEW
--------
This analysis searches for transposable elements in real bacterial genomes
by detecting inverted repeats (4-6 bp length) which are characteristic
markers of transposons.


BACTERIAL GENOMES ANALYZED
---------------------------
1. Escherichia coli (NC_000913.3) - ~4.6 Mbp
2. Bacillus subtilis (NC_000964.3) - ~4.2 Mbp
3. Mycoplasma genitalium (NC_000908.2) - ~0.58 Mbp


FILES
-----
  download_genomes.py                - Downloads genomes from NCBI
  find_transposons_real.py           - Finds inverted repeats in genomes
  bacterial_transposon_analysis.json - Results (output)


REQUIREMENTS
------------
  pip install biopython


USAGE
-----

Step 1: Download Bacterial Genomes

  python download_genomes.py

  This downloads 3 bacterial genomes from NCBI GenBank.


Step 2: Analyze Genomes for Transposons

  python find_transposons_real.py

  This:
  - Searches for inverted repeats (4-6 bp)
  - Identifies potential transposon insertion sites
  - Calculates transposon density per genome
  - Reports top 20 findings per genome


DETECTION ALGORITHM
-------------------

Inverted Repeat Detection:
  - Length range: 4-6 bp
  - Minimum spacing: 10 bp (between left and right TIR)
  - Maximum spacing: 100 bp
  - Filters out overlapping repeats

Example Transposon Structure:

  Sequence: ...ATGC[10-100bp]GCAT...
            ^^^^             ^^^^
           Left TIR      Right TIR (reverse complement)


WHY 4-6 BP?
-----------
Short inverted repeats (4-6 bp) are common in bacterial transposons,
particularly in insertion sequences (IS elements). While longer repeats
provide more specificity, shorter repeats are more sensitive for initial
screening of potential transposon sites.


INTERPRETATION
--------------
High density of inverted repeats may indicate:
  - Active transposon presence
  - Genome plasticity
  - Evolutionary adaptation mechanisms


