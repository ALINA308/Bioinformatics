# DNA Pattern Analysis - LAB13: BILCIURESCU ELENA-ALINA ( 1241 EA)

## About
This is my implementation for the LAB13 bioinformatics assignment. The program analyzes DNA sequences from gene promoters using a sliding window approach to compute C+G percentage and Kappa Index of Coincidence patterns.

## What I Did

I built a Python application that:
- Analyzes the full 82bp test sequence with a 30bp sliding window
- Calculates CG% (returns 29.27 for the test sequence)
- Calculates Kappa Index of Coincidence (returns 27.53 for the test sequence)
- Generates pattern charts showing how these values change across the sequence
- Analyzes multiple promoter sequences and compares their patterns
- Creates visualizations to compare different promoter patterns

## The Formulas I Used

**CG% (C+G percentage):**
```
CG% = (C_count + G_count) / sequence_length × 100
```
Returns 29.27 for the test sequence.

**IC (Kappa Index of Coincidence):**
```
IC = Σ(pi²) × 83
```
Where pi is the frequency of each nucleotide. Returns 27.53 for the test sequence.

**Center of Weight:**
```
Weighted average of positions using (CG% + IC) as weight
```

## How to Run

Just run the main script:
```bash
python ex1.py
```

This will:
- Analyze the 82bp test sequence
- Show that my functions return the correct values (29.27 and 27.53)
- Generate two charts:
  - `dna_pattern_analysis.png` - Shows the pattern for the test sequence
  - `dna_centers_of_weight.png` - Shows centers for all 5 promoters

## What Gets Generated

The program creates two PNG files:
1. **dna_pattern_analysis.png** - CG% and IC values across all 53 windows
2. **dna_centers_of_weight.png** - Comparison of 5 different promoter centers

## Results

For the test sequence:
- Full sequence CG%: 29.27 ✓
- Full sequence IC: 27.53 ✓
- 53 overlapping windows analyzed
- Center of weight at position 40.71bp

I also analyzed 4 additional promoters with different characteristics (AT-rich, GC-rich, balanced, A-dominant) and plotted all their centers on the second chart.

## Dependencies

Need to have these installed:
```bash
pip install matplotlib numpy
```

## Files in This Project

- `ex1.py` - Main program with all the functions
- `dna_pattern_analysis.png` - Pattern chart output
- `dna_centers_of_weight.png` - Centers comparison chart
- `README.md` - This file
