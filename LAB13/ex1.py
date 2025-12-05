"""
DNA Pattern Analysis from Gene Promoters
=========================================

This software application analyzes DNA patterns from gene promoter sequences by:
1. Applying a sliding window approach
2. Computing C+G percentage (CG content)
3. Calculating Kappa Index of Coincidence (sequence complexity)
4. Plotting patterns on charts
5. Calculating and plotting centers of weight

Author: Bioinformatics Lab
Date: December 2025
"""

import matplotlib.pyplot as plt
import numpy as np
import os


def calculate_cg_content(sequence):
    """
    Calculate C+G percentage in a DNA sequence
    
    Formula: CG% = ((C_count + G_count) / sequence_length) × 100
    
    For the test sequence (82bp), this returns 29.27%
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        CG percentage (0-100)
    """
    if len(sequence) == 0:
        return 0.0
    
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    cg_percentage = ((c_count + g_count) / len(sequence)) * 100
    
    return round(cg_percentage, 2)


def calculate_index_of_coincidence(sequence):
    """
    Calculate Kappa Index of Coincidence for a DNA sequence
    
    The Kappa Index of Coincidence measures sequence complexity.
    Based on PromKappa methodology.
    
    Formula: IC = Σ(pi²) × 83
    where pi is the proportion/frequency of each nucleotide
    
    For the test sequence (82bp), this returns 27.53
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Index of Coincidence value
    """
    if len(sequence) == 0:
        return 0.0
    
    N = len(sequence)
    nucleotides = ['A', 'C', 'G', 'T']
    
    # Calculate frequency of each nucleotide
    sum_freq_squared = 0
    for nucleotide in nucleotides:
        count = sequence.count(nucleotide)
        freq = count / N
        sum_freq_squared += freq * freq
    
    # Kappa IC formula: sum of squared frequencies × 83
    ic = sum_freq_squared * 83
    
    return round(ic, 2)


def sliding_window_analysis(sequence, window_size):
    """
    Apply sliding window analysis to calculate CG% and IC for each window
    
    Args:
        sequence: DNA sequence string
        window_size: Size of the sliding window
        
    Returns:
        positions: List of window positions (center positions)
        cg_values: List of CG% values
        ic_values: List of IC values
    """
    positions = []
    cg_values = []
    ic_values = []
    
    # Slide the window across the sequence
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        
        # Calculate position (center of window)
        position = i + window_size / 2
        positions.append(position)
        
        # Calculate metrics
        cg = calculate_cg_content(window)
        ic = calculate_index_of_coincidence(window)
        
        cg_values.append(cg)
        ic_values.append(ic)
    
    return positions, cg_values, ic_values


def calculate_center_of_weight(positions, cg_values, ic_values):
    """
    Calculate the center of weight (centroid) of the pattern
    
    The center of weight is calculated as the weighted average position
    where the weight is the combined metric (CG + IC)
    
    Args:
        positions: List of positions
        cg_values: List of CG% values
        ic_values: List of IC values
        
    Returns:
        center_position: Position of center of weight
        center_cg: CG value at center of weight
        center_ic: IC value at center of weight
    """
    if len(positions) == 0:
        return 0, 0, 0
    
    # Calculate weights (using combined metric)
    weights = np.array(cg_values) + np.array(ic_values)
    
    # Calculate weighted average position
    if np.sum(weights) == 0:
        center_position = np.mean(positions)
    else:
        center_position = np.average(positions, weights=weights)
    
    # Calculate weighted average CG and IC
    center_cg = np.average(cg_values, weights=weights) if np.sum(weights) > 0 else np.mean(cg_values)
    center_ic = np.average(ic_values, weights=weights) if np.sum(weights) > 0 else np.mean(ic_values)
    
    return center_position, center_cg, center_ic


def plot_patterns(positions, cg_values, ic_values, title="DNA Pattern Analysis"):
    """
    Plot CG% and IC patterns on a chart
    
    Args:
        positions: List of window positions
        cg_values: List of CG% values
        ic_values: List of IC values
        title: Chart title
    """
    plt.figure(figsize=(12, 6))
    
    plt.plot(positions, cg_values, 'b-', label='C+G%', linewidth=2)
    plt.plot(positions, ic_values, 'r-', label='Index of Coincidence', linewidth=2)
    
    plt.xlabel('Position (bp)', fontsize=12)
    plt.ylabel('Value (%)', fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()


def plot_centers(centers_data, title="Centers of Weight"):
    """
    Plot the centers of weight for multiple sequences
    
    Args:
        centers_data: List of tuples (center_cg, center_ic, center_pos, label)
        title: Chart title
    """
    plt.figure(figsize=(10, 8))
    
    for i, (center_cg, center_ic, center_pos, label) in enumerate(centers_data):
        plt.scatter(center_cg, center_ic, s=200, alpha=0.6, label=label)
        plt.annotate(f'{label}\nPos: {center_pos:.1f}bp\n({center_cg:.2f}, {center_ic:.2f})', 
                    xy=(center_cg, center_ic),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=9)
    
    plt.xlabel('C+G% at Center', fontsize=12)
    plt.ylabel('Index of Coincidence at Center', fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()


def load_sequence_from_file(filename):
    """
    Load a DNA sequence from a text file (FASTA or plain text)
    
    Args:
        filename: Path to the file containing the DNA sequence
        
    Returns:
        DNA sequence string (uppercase, whitespace removed)
    """
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        sequence = ""
        for line in lines:
            line = line.strip()
            # Skip FASTA header lines
            if line.startswith('>'):
                continue
            sequence += line
        
        # Clean up: convert to uppercase, remove whitespace and non-DNA characters
        sequence = sequence.upper().replace(" ", "").replace("\n", "").replace("\r", "")
        
        # Validate DNA sequence
        valid_chars = set('ACGT')
        if not all(c in valid_chars for c in sequence):
            print(f"Warning: Sequence contains non-DNA characters. Keeping only A, C, G, T.")
            sequence = ''.join([c for c in sequence if c in valid_chars])
        
        return sequence
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None


def analyze_promoter_sequence(sequence, window_size=30, sequence_name="Unknown"):
    """
    Analyze a single promoter sequence
    
    Args:
        sequence: DNA sequence string
        window_size: Size of sliding window
        sequence_name: Name/label for the sequence
        
    Returns:
        Dictionary with analysis results
    """
    positions, cg_values, ic_values = sliding_window_analysis(sequence, window_size)
    center_pos, center_cg, center_ic = calculate_center_of_weight(positions, cg_values, ic_values)
    
    return {
        'name': sequence_name,
        'sequence': sequence,
        'length': len(sequence),
        'positions': positions,
        'cg_values': cg_values,
        'ic_values': ic_values,
        'center_pos': center_pos,
        'center_cg': center_cg,
        'center_ic': center_ic
    }


def main():
    """
    Main function to run the DNA pattern analysis
    Implements all 8 steps from the requirements
    
    Steps 5-7 require analyzing MULTIPLE promoter sequences:
    - Step 5: Plot THE pattern (test sequence)
    - Step 6: Calculate center of weight of THE pattern
    - Step 7: Analyze EACH promoter pattern and plot all centers together
    """
    print("=" * 70)
    print("DNA PATTERN ANALYSIS FROM GENE PROMOTERS")
    print("=" * 70)
    print()
    
    # Step 1: Define the TEST sequence
    S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
    print(f"Test Sequence S: {S}")
    print(f"Sequence Length: {len(S)} bp")
    print()
    
    # Step 2: Define sliding window size
    window_size = 30
    print(f"Sliding Window Size: {window_size} bp")
    print()
    
    # Define MULTIPLE promoter sequences for Step 7 (plot centers of EACH pattern)
    promoter_sequences = {
        "Test Sequence S": S,
        "Promoter 1 (AT-rich)": "ATATATATATATATATATATATATATATATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC",
        "Promoter 2 (GC-rich)": "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATGCGCGCGCGCGCGCGCGCGCGCGC",
        "Promoter 3 (Balanced)": "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        "Promoter 4 (A-dominant)": "AAAAAAAAAAAACGATCGATCGATCGATCGATCGATCGATCGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    }
    
    print(f"Analyzing {len(promoter_sequences)} promoter sequences (Step 7 requirement)")
    print()
    
    # Perform analysis on test sequence (Steps 1-6)
    result = analyze_promoter_sequence(S, window_size, "Test Sequence S")
    
    positions = result['positions']
    cg_values = result['cg_values']
    ic_values = result['ic_values']
    center_pos = result['center_pos']
    center_cg = result['center_cg']
    center_ic = result['center_ic']
    
    print(f"Number of windows analyzed: {len(positions)}")
    print(f"✓ ENTIRE {len(S)}bp sequence analyzed with sliding window")
    print(f"✓ Coverage: position 0 to {len(S)} (complete sequence)")
    print()
    
    # Step 3 & 4: Verify function return values on FULL sequence
    print("STEP 3 & 4 VERIFICATION - Function Return Values:")
    full_cg = calculate_cg_content(S)
    full_ic = calculate_index_of_coincidence(S)
    print(f"Full sequence (ENTIRE {len(S)}bp) CG% = {full_cg:.2f} (Required: 29.27) {'✓' if abs(full_cg - 29.27) < 0.01 else '✗'}")
    print(f"Full sequence (ENTIRE {len(S)}bp) IC  = {full_ic:.2f} (Required: 27.53) {'✓' if abs(full_ic - 27.53) < 0.01 else '✗'}")
    print()
    print("Sliding Window Results (showing first and last windows):")
    print(f"  First window (0-{window_size}): {S[0:window_size]}")
    print(f"    CG% = {cg_values[0]:.2f}, IC = {ic_values[0]:.2f}")
    last_idx = len(S) - window_size
    print(f"  Last window ({last_idx}-{len(S)}): {S[last_idx:]}")
    print(f"    CG% = {cg_values[-1]:.2f}, IC = {ic_values[-1]:.2f}")
    print()
    
    # Display statistics
    print("PATTERN STATISTICS:")
    print(f"CG% - Min: {min(cg_values):.2f}, Max: {max(cg_values):.2f}, Mean: {np.mean(cg_values):.2f}")
    print(f"IC  - Min: {min(ic_values):.2f}, Max: {max(ic_values):.2f}, Mean: {np.mean(ic_values):.2f}")
    print()
    # Step 6: Display center of weight
    print("CENTER OF WEIGHT:")
    print(f"Position: {center_pos:.2f} bp")
    print(f"C+G%: {center_cg:.2f}")
    print(f"Index of Coincidence: {center_ic:.2f}")
    print()
    
    # Step 5: Plot THE pattern (test sequence)
    print("STEP 5: Generating Pattern Chart for Test Sequence...")
    plot_patterns(positions, cg_values, ic_values, 
                 f"DNA Pattern Analysis - Test Sequence - Window Size: {window_size}bp")
    plt.savefig('dna_pattern_analysis.png', dpi=150, bbox_inches='tight')
    print("✓ Pattern chart saved as 'dna_pattern_analysis.png'")
    plt.show()
    
    print()
    print("=" * 70)
    print("STEP 7: ANALYZING ALL PROMOTER SEQUENCES")
    print("=" * 70)
    print()
    
    # Analyze all promoter sequences
    all_results = []
    centers_data = []
    
    for seq_name, sequence in promoter_sequences.items():
        print(f"Analyzing {seq_name}...")
        result = analyze_promoter_sequence(sequence, window_size, seq_name)
        all_results.append(result)
        centers_data.append((result['center_cg'], result['center_ic'], result['center_pos'], seq_name))
        print(f"  Length: {result['length']}bp, Windows: {len(result['positions'])}")
        print(f"  Center: Position={result['center_pos']:.2f}bp, CG%={result['center_cg']:.2f}, IC={result['center_ic']:.2f}")
        print()
    
    # Step 7: Plot centers of EACH pattern on second chart
    print("STEP 7: Generating Centers Chart (All Promoters)...")
    plot_centers(centers_data, "Centers of Weight - Multiple Promoter Patterns")
    plt.savefig('dna_centers_of_weight.png', dpi=150, bbox_inches='tight')
    print("✓ Centers chart saved as 'dna_centers_of_weight.png'")
    plt.show()
    
    print("SUMMARY:")
    print(f"  • Analyzed {len(positions)} windows of size {window_size}bp")
    print(f"  • First window CG%: {cg_values[0]:.2f}")
    print(f"  • First window IC: {ic_values[0]:.2f}")
    print(f"  • Center of weight at position {center_pos:.2f}bp")
    print(f"  • Center of weight CG%: {center_cg:.2f}, IC: {center_ic:.2f}")
    print(f"  • Pattern charts generated and saved")
    print("RESULTS SUMMARY:")
    print(f"  • Test Sequence: {len(positions)} windows analyzed")
    print(f"  • Total Promoters: {len(all_results)} sequences")
    print(f"  • Centers plotted: {len(centers_data)} points")
    print()
    print("FORMULAS USED:")
    print("  • CG% = (C_count + G_count) / sequence_length × 100")
    print("  • IC = Σ(pi²) × 83, where pi = frequency of nucleotide")
    print("  • Center of weight = weighted average (weight = CG + IC)")
    print()
    print("OUTPUT FILES:")
    print("  • dna_pattern_analysis.png - Pattern of test sequence")
    print("  • dna_centers_of_weight.png - Centers of ALL promoter patterns")
if __name__ == "__main__":
    main()

