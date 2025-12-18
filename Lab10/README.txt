================================================================================
LAB 10 - DNA MOTIF FINDING FOR EXON-INTRON BORDER DETECTION
================================================================================

AUTHOR: Bilciurescu Elena-Alina 1241 EA
DATE: December 18, 2025

================================================================================
OVERVIEW
================================================================================

This lab implements a complete DNA motif finding application for identifying
exon-intron boundaries in DNA sequences. The application uses log-likelihood
scoring to detect splice site signals in both test sequences and real
influenza virus genomes.

================================================================================
FILES INCLUDED
================================================================================

1. complete_solution.py
   - Complete implementation of exercises 1-5
   - Builds all matrices (Count, Weight, Frequencies, Log-likelihoods)
   - Analyzes test sequence S for exon-intron borders
   - Shows step-by-step calculations

2. scan_influenza_genomes.py
   - Adapted application for scanning influenza genomes
   - Scans 10 influenza virus genomes for motif patterns
   - Generates visualization charts for each genome
   - Identifies most likely functional motif locations

3. influenza_genomes.fasta
   - FASTA file containing 10 influenza virus genome sequences
   - Downloaded from NCBI database
   - Total sequences: 10 segments from H3N2 and H1N1 strains

4. genome_*_motif_chart.png (10 files)
   - Visualization charts showing motif analysis for each genome
   - Upper panel: Log-likelihood score distribution across genome
   - Lower panel: Genome map with marked motif locations
   - Color-coded by score intensity (red=high, orange=medium, yellow=low)

================================================================================
EXERCISE 1-5: MOTIF FINDING BASICS
================================================================================

TASK:
Analyze sequence S using the exon-intron boundary motif from 9 known sequences.

APPROACH:
1. Count Matrix: Count nucleotides at each position
2. Weight Matrix: Add pseudocounts to avoid zero probabilities
3. Frequency Matrix: Convert counts to probabilities
4. Log-likelihood Matrix: Compare to null model (random DNA)
5. Sliding Window: Score all possible 9-nucleotide windows in sequence S

RESULTS:
- Test sequence: CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA (38 bp)
- Windows analyzed: 30
- Significant motifs found: 5 (with positive scores)
- Best match: 'AACGTAATC' at position 9 (score: 3.516)
- CONCLUSION: YES - Sequence S contains exon-intron border signal!

TO RUN:
python complete_solution.py

================================================================================
INFLUENZA GENOME ANALYSIS
================================================================================

TASK:
Download 10 influenza genomes and scan for exon-intron boundary motifs.
Create charts showing most likely locations of functional motifs.

GENOMES ANALYZED:
1. CY121496.1 - Influenza A (H3N2) HA gene (1,732 bp) - 259 motifs
2. CY121497.1 - Influenza A (H3N2) M1/M2 genes (1,003 bp) - 157 motifs
3. CY121498.1 - Influenza A (H3N2) NA gene (1,443 bp) - 241 motifs
4. CY121499.1 - Influenza A (H3N2) segment 4 (1,541 bp) - 232 motifs
5. CY121500.1 - Influenza A (H3N2) segment 5 (865 bp) - 129 motifs
6. CY121501.1 - Influenza A (H3N2) segment 6 (2,208 bp) - 336 motifs
7. CY121502.1 - Influenza A (H3N2) segment 7 (2,316 bp) - 329 motifs
8. CY121503.1 - Influenza A (H3N2) segment 8 (2,309 bp) - 373 motifs
9. KF848938.1 - Influenza A (H1N1) segment 1 (1,665 bp) - 248 motifs
10. KF848939.1 - Influenza A (H1N1) segment 2 (586 bp) - 81 motifs

KEY FINDINGS:
- Total windows analyzed: 15,638
- Total significant motifs (score > 0): 2,385
- Highest scoring motif: 5.384 in genome CY121503.1 (segment 8)
- Average motifs per genome: 238.5

CHARTS GENERATED:
Each chart shows:
- Score distribution plot along entire genome
- Highlighted regions with positive scores (potential splice sites)
- Top 10 strongest motifs marked with red stars
- Genome map with color-coded motif intensity
- Statistical summary (total motifs, max score, etc.)

TO RUN:
python scan_influenza_genomes.py

NOTE: Requires matplotlib library:
pip install matplotlib

================================================================================
ALGORITHM DETAILS
================================================================================

LOG-LIKELIHOOD SCORING:
For each 9-nucleotide window:
  Score = Sum of log(P(nucleotide_i at position_i) / 0.25)

  where:
  - P(nucleotide_i at position_i) comes from frequency matrix
  - 0.25 is the null model (random DNA probability)

INTERPRETATION:
  - Positive score: Window is MORE similar to motif than random
  - Negative score: Window is LESS similar to motif than random
  - Score > 2: Strong candidate for functional motif
  - Score > 3: Very strong candidate

VISUALIZATION:
  - Dark red: High confidence (score > 2)
  - Red: Medium confidence (1 < score <= 2)
  - Orange: Low confidence (0 < score <= 1)
  - Gray: No signal (score <= 0)

================================================================================
BIOLOGICAL SIGNIFICANCE
================================================================================

EXON-INTRON BOUNDARIES:
The motif pattern detected represents splice sites - critical sequences that
mark the boundaries between coding (exons) and non-coding (introns) regions
in pre-mRNA.

CONSENSUS SEQUENCE:
From the 9 training sequences, the conserved pattern is:
  Position:  1  2  3  4  5  6  7  8  9
  Consensus: N  A  G  G  T  A  A  N  N

Where position 4 (G) and position 5 (T) are highly conserved (GT dinucleotide
is the canonical 5' splice site donor sequence).

INFLUENZA VIRUS:
Influenza has a segmented RNA genome. The presence of multiple high-scoring
motifs in each segment suggests active splicing mechanisms, which is important
for viral protein expression and regulation.

================================================================================
TECHNICAL NOTES
================================================================================

PSEUDOCOUNTS:
Added 1 pseudocount to each cell in count matrix to avoid log(0) errors.
This Laplace smoothing ensures all probabilities are non-zero.

SEQUENCE HANDLING:
Ambiguous nucleotides (N, R, Y, etc.) are assigned very low scores (-10)
to prevent false positives from uncertain sequence data.

PERFORMANCE:
- Average scan time: ~0.5 seconds per genome
- Memory efficient: Processes sequences sequentially
- Handles genomes from 586 bp to 2,316 bp

STATISTICAL THRESHOLD:
Using score > 0 as cutoff ensures we only report windows more likely than
random chance. For high-confidence predictions, use score > 2.

================================================================================
RESULTS SUMMARY
================================================================================

EXERCISE 5 ANSWER:
YES - The test sequence S contains signals indicating an exon-intron border.
The motif 'AACGTAATC' at position 9 has a log-likelihood score of 3.516,
which is 33.6 times more likely than a random sequence.

INFLUENZA ANALYSIS:
All 10 influenza genomes show multiple potential splice sites, with an
average of 238.5 significant motifs per genome. The highest scoring motif
(5.384) was found in segment 8 of the H3N2 strain at position 2102.


