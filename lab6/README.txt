Author: Bilciurescu Elena-Alina

PROJECT: Gel Electrophoresis Simulation

DESCRIPTION:
This project simulates DNA gel electrophoresis, a fundamental analysis method
used in life sciences to determine the relative sizes of DNA fragments.

IMPLEMENTATION:
1. DNA Sequence Input: Reads a DNA sequence from NCBI (1000-3000 nucleotides)
   stored in FASTA format.

2. Random Sampling: Extracts 10 random DNA fragments from the original sequence,
   with fragment lengths ranging from 100 to 3000 base pairs.

3. Data Storage: Stores all extracted samples in an array structure with their
   properties (ID, sequence, length, start position).

4. Electrophoresis Simulation: Simulates DNA fragment migration on an
   electrophoresis gel based on fragment length (molecular weight proxy).
   - Short fragments = less friction = faster migration = travel further
   - Long fragments = more friction = slower migration = travel less distance

5. Visualization: Generates a visual representation of the gel with:
   - DNA ladder (reference markers at 500, 1500, 3000 bp)
   - Sample lanes showing band positions
   - Migration distances inversely proportional to fragment length

OUTPUT:
- Console output with fragment analysis table
- PNG image of simulated gel electrophoresis (gel_electrophoresis.png)

HOW TO RUN:
python gel_electrophoresis.py
