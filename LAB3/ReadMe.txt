Author:
     Bilciurescu Elena-Alina
(Worked alone on this project)

Project: DNA Melting Temperature Calculator (Lab 3)

================================================================================
EXERCISE 1: DNA MELTING TEMPERATURE CALCULATOR
================================================================================

Description:
This application calculates the melting temperature (Tm) of a DNA sequence using two different formulas:

1. Simple Formula: Tm = 4(G + C) + 2(A + T) degrees Celsius
2. Salt-Adjusted Formula: Tm = 81.5 + 16.6(log10([Na+])) + 0.41*(%GC) - 600/length

How to Run:

Step 1: Open a terminal or command prompt

Step 2: Navigate to the directory containing ex1.py
        Example: cd c:\Users\alina\Desktop\LAB3

Step 3: Run the program using Python
        Command: python ex1.py

Step 4: Follow the on-screen prompts:
        a) Enter your DNA sequence (e.g., ATCGATCG, GCGCGCGC)
        b) Choose whether to use default Na+ concentration (50 mM) by entering 'y' or 'n'
        c) If you chose 'n', enter your desired Na+ concentration in mM

Step 5: The program will display:
        - Your DNA sequence and its length
        - Melting temperature using the simple formula
        - Melting temperature using the salt-adjusted formula

Input:
- DNA sequence (string containing A, T, G, C nucleotides - case insensitive)
- Optional: Na+ concentration in mM for salt-adjusted formula (default: 50 mM)

Output:
- Melting temperature in Celsius using both methods

Features:
- Input validation for DNA sequences (only A, T, G, C allowed)
- Both calculation methods implemented
- Customizable Na+ concentration for salt-adjusted formula
- User-friendly interface with formatted output
- Handles both uppercase and lowercase DNA sequences


================================================================================
EXERCISE 2: SLIDING WINDOW MELTING TEMPERATURE ANALYZER
================================================================================

Description:
This application uses the sliding window method to calculate melting temperature
across a DNA sequence. It slides an 8-base pair window along the sequence and
calculates the Tm for each position, allowing you to identify regions of varying
thermal stability. The results are visualized in a chart/plot.

Melting Temperature Formula (Wallace Rule for sequences ≤14 bp):
- Tm = 2(A + T) + 4(G + C) degrees Celsius

For longer sequences:
- Tm = 64.9 + 41*(G+C-16.4)/(A+T+G+C)

How to Run:

Method 1 - Interactive Mode:
Step 1: Open a terminal or command prompt
Step 2: Navigate to the directory containing ex2.py
        Example: cd c:\Users\alina\Desktop\LAB3
Step 3: Run the program: python ex2.py
Step 4: Enter the path to your FASTA file when prompted
        Example: sample_sequence.fasta
Step 5: View the console results and the generated chart

Method 2 - Command-Line Mode:
Step 1: Open a terminal or command prompt
Step 2: Navigate to the directory
Step 3: Run with FASTA file as argument:
        Command: python ex2.py sample_sequence.fasta

Prerequisites:
- Python 3.x
- matplotlib library (install with: pip install matplotlib)

Input:
- FASTA file containing a DNA sequence
- Sample file provided: sample_sequence.fasta

Output:
- Console output: Table showing position, window sequence, and melting temperature for each window
- Console output: Statistics (Average, Minimum, and Maximum Tm values)
- Chart/Plot: Visual representation of the melting temperature profile across the sequence

Features:
- FASTA file parser
- Sliding window analysis (window size: 8 bp)
- Melting temperature calculation for each window position
- Comprehensive statistics (average, min, max Tm)
- Graphical visualization with matplotlib
- Blue line showing Tm values at each position
- Red dashed line indicating average Tm
- Grid and proper axis labels for easy interpretation
- Handles sequences of any length
- Command-line and interactive mode support

Example Output (Console):
Position    Window Sequence    Tm (C)
--------------------------------------
0          AAAAAAAA           16.00
1          AAAAAATT           16.00
2          AAAAATTT           16.00
...

STATISTICS:
  Average Tm: 24.50 C
  Minimum Tm: 16.00 C
  Maximum Tm: 32.00 C

Chart:
A matplotlib window will open displaying a line graph showing the melting
temperature profile across the entire sequence, with peaks in GC-rich regions
and valleys in AT-rich regions.
