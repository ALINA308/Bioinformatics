Bioinformatics Lab 4 - Codon Frequency Comparison

Author
Alina Bilciurescu

Project Description
This project solves a bioinformatics problem by downloading FASTA files for the COVID-19 and Influenza genomes from NCBI, analyzing codon frequencies, creating comparative charts, and identifying top amino acids. The analysis includes:

- Downloading genome sequences from NCBI.
- Counting codon frequencies across the entire sequences.
- Generating charts for top 10 codons in each virus and a comparison chart.
- Translating sequences to amino acids and counting frequencies.
- Outputting top 3 amino acids for each genome in the console.
- Formulating prompts for an AI to identify foods lacking the top amino acids.

Files in the Project
- ex1.py: Script for translating DNA/RNA sequences to amino acids using the standard genetic code.
- ex2.py: Main script that downloads FASTA files, processes sequences, counts codons and amino acids, generates charts, and prints results.
- covid19.fasta: Downloaded FASTA file for SARS-CoV-2 (COVID-19) genome.
- influenza.fasta: Downloaded FASTA file for Influenza A virus genome.
- covid_codons.png: Chart of top 10 codons in COVID-19.
- flu_codons.png: Chart of top 10 codons in Influenza.
- comparison_codons.png: Comparison chart of top codons between COVID-19 and Influenza.


Problems Solved
1. Download FASTA Files: Retrieve genome sequences from NCBI using Entrez efetch API.
2. Sequence Processing: Normalize sequences (convert T to U for RNA), parse FASTA format.
3. Codon Counting: Count frequencies of all codons in the sequences.
4. Amino Acid Analysis: Translate codons to amino acids and count frequencies.
5. Visualization: Create bar charts for top codons and comparisons using Matplotlib.
6. Output Results: Print top 10 codons and top 3 amino acids to console.
7. AI Prompt Formulation: Generate questions about foods lacking specific amino acids.

How to Run
1. Ensure Python 3.x is installed with required libraries: urllib, collections, matplotlib.
2. Run python ex2.py to execute the analysis. It will download files, process data, generate charts, and print results.
3. View generated PNG files for charts.

Dependencies
- Python 3.x
- urllib.request (built-in)
- collections.Counter (built-in)
- matplotlib


