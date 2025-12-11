"""
Exercise 2: Local Alignment of Influenza and COVID-19 Genomes
Uses Smith-Waterman algorithm with chunking for large sequences
"""

import json


def smith_waterman(seq1, seq2, match_score=2, mismatch_penalty=-1, gap_penalty=-2):
    n = len(seq1)
    m = len(seq2)

    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i-1] == seq2[j-1]:
                diagonal_score = score_matrix[i-1][j-1] + match_score
            else:
                diagonal_score = score_matrix[i-1][j-1] + mismatch_penalty

            up_score = score_matrix[i-1][j] + gap_penalty
            left_score = score_matrix[i][j-1] + gap_penalty

            score_matrix[i][j] = max(0, diagonal_score, up_score, left_score)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    return max_score, max_pos, score_matrix


def traceback_smith_waterman(seq1, seq2, score_matrix, max_pos, match_score=2, mismatch_penalty=-1, gap_penalty=-2):
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = max_pos

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        current_score = score_matrix[i][j]

        if seq1[i-1] == seq2[j-1]:
            diagonal_score = score_matrix[i-1][j-1] + match_score
        else:
            diagonal_score = score_matrix[i-1][j-1] + mismatch_penalty

        if current_score == diagonal_score:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif current_score == score_matrix[i-1][j] + gap_penalty:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        elif current_score == score_matrix[i][j-1] + gap_penalty:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        else:
            break

    start_pos = (i, j)
    return aligned_seq1, aligned_seq2, start_pos


def chunked_local_alignment(seq1, seq2, chunk_size=500, overlap=100, match_score=2, mismatch_penalty=-1, gap_penalty=-2):
    results = []
    num_chunks = (len(seq1) - overlap) // (chunk_size - overlap) + 1

    print(f"\nDividing sequence 1 ({len(seq1)} bp) into {num_chunks} chunks of {chunk_size} bp with {overlap} bp overlap")
    print(f"Aligning against sequence 2 ({len(seq2)} bp)\n")

    for chunk_idx in range(num_chunks):
        start_idx = chunk_idx * (chunk_size - overlap)
        end_idx = min(start_idx + chunk_size, len(seq1))

        chunk = seq1[start_idx:end_idx]

        if len(chunk) < 50:
            continue

        print(f"Processing chunk {chunk_idx + 1}/{num_chunks} (positions {start_idx}-{end_idx})...", end=" ")

        max_score, max_pos, score_matrix = smith_waterman(chunk, seq2, match_score, mismatch_penalty, gap_penalty)

        aligned_seq1, aligned_seq2, start_pos = traceback_smith_waterman(
            chunk, seq2, score_matrix, max_pos, match_score, mismatch_penalty, gap_penalty
        )

        matches = sum(1 for k in range(len(aligned_seq1)) if aligned_seq1[k] == aligned_seq2[k])
        identity = (matches / len(aligned_seq1) * 100) if len(aligned_seq1) > 0 else 0

        result = {
            'chunk_index': chunk_idx,
            'seq1_start': start_idx,
            'seq1_end': end_idx,
            'seq2_align_start': start_pos[1],
            'seq2_align_end': max_pos[1],
            'alignment_score': max_score,
            'alignment_length': len(aligned_seq1),
            'matches': matches,
            'identity_percentage': identity,
            'aligned_seq1': aligned_seq1,
            'aligned_seq2': aligned_seq2
        }

        results.append(result)
        print(f"Score: {max_score}, Length: {len(aligned_seq1)}, Identity: {identity:.2f}%")

    return results


def visualize_genome_similarities(results, seq1_length, seq2_length):
    print("\n" + "=" * 100)
    print("GENOME SIMILARITY VISUALIZATION")
    print("=" * 100)

    display_width = 100

    print(f"\nSequence 1 (Query): {seq1_length} bp")
    print(f"Sequence 2 (Reference): {seq2_length} bp")
    print(f"\nAlignment regions (scaled to {display_width} characters):\n")

    visualization = [' '] * display_width

    for result in results:
        if result['alignment_score'] > 0:
            start_display = int((result['seq1_start'] / seq1_length) * display_width)
            end_display = int((result['seq1_end'] / seq1_length) * display_width)

            if result['identity_percentage'] >= 70:
                marker = '#'
            elif result['identity_percentage'] >= 50:
                marker = '+'
            else:
                marker = '.'

            for pos in range(start_display, min(end_display + 1, display_width)):
                visualization[pos] = marker

    print("Query:    [", end="")
    print(''.join(visualization), end="")
    print("]")

    print("\nLegend: # = High similarity (>=70%), + = Medium similarity (>=50%), . = Low similarity (<50%)\n")

    print("=" * 100)
    print("DETAILED ALIGNMENT STATISTICS")
    print("=" * 100)
    print(f"{'Chunk':<8}{'Query Pos':<20}{'Ref Pos':<20}{'Score':<10}{'Length':<10}{'Identity':<12}")
    print("-" * 100)

    total_aligned_bases = 0
    total_matches = 0

    for result in results:
        if result['alignment_score'] > 0:
            print(f"{result['chunk_index']+1:<8}"
                  f"{result['seq1_start']}-{result['seq1_end']:<18}"
                  f"{result['seq2_align_start']}-{result['seq2_align_end']:<18}"
                  f"{result['alignment_score']:<10}"
                  f"{result['alignment_length']:<10}"
                  f"{result['identity_percentage']:.2f}%")
            total_aligned_bases += result['alignment_length']
            total_matches += result['matches']

    print("-" * 100)
    print(f"\nTotal aligned bases: {total_aligned_bases}")
    print(f"Total matches: {total_matches}")
    if total_aligned_bases > 0:
        print(f"Overall identity: {total_matches / total_aligned_bases * 100:.2f}%")


def save_alignment_results(results, filename):
    json_results = []
    for result in results:
        json_result = {
            'chunk_index': result['chunk_index'],
            'seq1_start': result['seq1_start'],
            'seq1_end': result['seq1_end'],
            'seq2_align_start': result['seq2_align_start'],
            'seq2_align_end': result['seq2_align_end'],
            'alignment_score': result['alignment_score'],
            'alignment_length': result['alignment_length'],
            'matches': result['matches'],
            'identity_percentage': result['identity_percentage'],
            'aligned_seq1_preview': result['aligned_seq1'][:100] + "..." if len(result['aligned_seq1']) > 100 else result['aligned_seq1'],
            'aligned_seq2_preview': result['aligned_seq2'][:100] + "..." if len(result['aligned_seq2']) > 100 else result['aligned_seq2']
        }
        json_results.append(json_result)

    with open(filename, 'w') as f:
        json.dump(json_results, f, indent=2)

    print(f"\nAlignment results saved to {filename}")


def parse_fasta(filename):
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                sequence += line.upper()
    return sequence


if __name__ == "__main__":
    print("=" * 100)
    print("LOCAL ALIGNMENT OF INFLUENZA AND COVID-19 GENOMES")
    print("=" * 100)

    print("\nLoading genomes...")
    influenza_seq = parse_fasta("influenza_genome.fasta")
    covid_seq = parse_fasta("covid19_genome.fasta")

    print(f"Influenza genome: {len(influenza_seq)} bp")
    print(f"COVID-19 genome: {len(covid_seq)} bp")

    results = chunked_local_alignment(
        seq1=influenza_seq,
        seq2=covid_seq,
        chunk_size=500,
        overlap=100,
        match_score=2,
        mismatch_penalty=-1,
        gap_penalty=-2
    )

    visualize_genome_similarities(results, len(influenza_seq), len(covid_seq))
    save_alignment_results(results, "alignment_results.json")

    print("\n" + "=" * 100)
    print("EXAMPLE ALIGNMENTS (Top 3 by score)")
    print("=" * 100)

    sorted_results = sorted(results, key=lambda x: x['alignment_score'], reverse=True)

    for idx, result in enumerate(sorted_results[:3], 1):
        print(f"\n--- Example {idx} ---")
        print(f"Chunk {result['chunk_index']+1}: Query positions {result['seq1_start']}-{result['seq1_end']}")
        print(f"Reference positions: {result['seq2_align_start']}-{result['seq2_align_end']}")
        print(f"Score: {result['alignment_score']}, Identity: {result['identity_percentage']:.2f}%")

        seq1_display = result['aligned_seq1'][:100]
        seq2_display = result['aligned_seq2'][:100]

        print(f"\nQuery: {seq1_display}")

        match_str = ""
        for i in range(len(seq1_display)):
            if seq1_display[i] == seq2_display[i]:
                match_str += "|"
            elif seq1_display[i] == "-" or seq2_display[i] == "-":
                match_str += " "
            else:
                match_str += "x"
        print(f"       {match_str}")
        print(f"Ref:   {seq2_display}")

        if len(result['aligned_seq1']) > 100:
            print(f"... (showing first 100 of {len(result['aligned_seq1'])} aligned bases)")

    print("\n" + "=" * 100)
