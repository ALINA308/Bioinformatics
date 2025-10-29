================================================================================
DNA SEQUENCE RECONSTRUCTION FROM RANDOM SAMPLES - LABORATORY EXERCISE 5
================================================================================

AUTHOR(S):
----------
Bilciurescu Elena Alina 1241EA


EXERCISE DESCRIPTION:
---------------------
This project solves a bioinformatics problem involving DNA sequence
reconstruction from random fragments:

1. Take an arbitrary DNA sequence from NCBI (National Center for
   Biotechnology Information), between 1000 and 3000 nucleotides.

2. Take 2000 random samples from this sequence, of about 100-150 bases each.

3. Store these samples in an array.

4. Rebuild the original DNA sequence using these random samples.


IMPLEMENTATION:
---------------
File: ex1.py

The solution uses:
- A real human insulin DNA sequence (2000 bp) from NCBI
- Random sampling to generate 2000 fragments (100-150 bp each)
- Greedy overlap assembly algorithm for reconstruction
- Accuracy metrics to evaluate reconstruction quality


ALGORITHM APPROACH:
-------------------
The reconstruction uses a greedy overlap-based assembly algorithm:
1. Remove duplicate fragments
2. Start with the longest fragment
3. Iteratively find the best overlapping fragment (minimum 20 bp overlap)
4. Extend the sequence by merging fragments
5. Continue until no more overlaps are found


MAIN PROBLEMS WITH THIS ALGORITHM APPROACH:
--------------------------------------------

1. REPEAT SEQUENCES (Most Critical Problem):
   - DNA sequences often contain repetitive regions
   - The greedy algorithm cannot distinguish between different copies of
     the same repeat
   - This leads to misassembly, collapsed repeats, or incorrect connections
   - Example: When encountering "ATGATGATGATG", the algorithm cannot
     determine which "ATG" copy connects to which fragment

2. AMBIGUOUS OVERLAPS:
   - Multiple fragments may have similar overlap scores
   - Greedy algorithm picks the first/best option without considering
     global optimization
   - No backtracking mechanism if a wrong choice is made early
   - Local optimal choices don't guarantee global optimal solution

3. COVERAGE GAPS:
   - Random sampling may miss some genomic regions entirely
   - Even with 2000 samples, statistical gaps are possible
   - Uncovered regions lead to incomplete or fragmented reconstruction

4. COMPUTATIONAL COMPLEXITY:
   - Greedy overlap approach has O(n²) or worse time complexity for n fragments
   - Not scalable for large genomes or whole-genome assembly
   - Modern assemblers use more efficient data structures (e.g., De Bruijn
     graphs)

5. NO ERROR CORRECTION:
   - Algorithm doesn't handle sequencing errors or mutations
   - Single base pair difference can break overlap detection
   - Real sequencing data typically has ~1% error rate

6. ORIENTATION ISSUES:
   - This implementation only checks forward strand direction
   - Real DNA sequencing produces fragments from both strands
   - Should check reverse complement orientations as well


BETTER APPROACHES:
------------------
Modern genome assembly uses more sophisticated methods:
- De Bruijn graph assembly (Velvet, SPAdes assemblers)
- Overlap-Layout-Consensus (OLC) with graph optimization
- Paired-end reads for scaffolding and gap resolution
- Multiple sequence alignment algorithms
- Error correction and consensus calling
- Hybrid approaches combining short and long reads


RESULTS:
--------
The algorithm can successfully reconstruct sequences without repeats when
coverage is adequate, but struggles with the fundamental challenges listed
above. The accuracy metrics show the quality of reconstruction.


FILES:
------
- ex1.py        : Main implementation
- readme.txt    : This file


REQUIREMENTS:
-------------
- Python 3.x
- random module (standard library)


HOW TO RUN:
-----------
python ex1.py


================================================================================
