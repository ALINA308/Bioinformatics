"""
Exercise 1: Needleman-Wunsch Global Alignment Algorithm
"""

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-2):
    n = len(seq1)
    m = len(seq2)

    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    for i in range(n + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(m + 1):
        score_matrix[0][j] = j * gap_penalty

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i-1] == seq2[j-1]:
                diagonal_score = score_matrix[i-1][j-1] + match_score
            else:
                diagonal_score = score_matrix[i-1][j-1] + mismatch_penalty

            up_score = score_matrix[i-1][j] + gap_penalty
            left_score = score_matrix[i][j-1] + gap_penalty

            score_matrix[i][j] = max(diagonal_score, up_score, left_score)

    aligned_seq1 = ""
    aligned_seq2 = ""
    i = n
    j = m

    while i > 0 or j > 0:
        current_score = score_matrix[i][j]

        if i > 0 and j > 0:
            if seq1[i-1] == seq2[j-1]:
                diagonal_score = score_matrix[i-1][j-1] + match_score
            else:
                diagonal_score = score_matrix[i-1][j-1] + mismatch_penalty

            if current_score == diagonal_score:
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                i -= 1
                j -= 1
                continue

        if i > 0:
            up_score = score_matrix[i-1][j] + gap_penalty
            if current_score == up_score:
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                i -= 1
                continue

        if j > 0:
            left_score = score_matrix[i][j-1] + gap_penalty
            if current_score == left_score:
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                j -= 1
                continue

    final_score = score_matrix[n][m]
    return aligned_seq1, aligned_seq2, score_matrix, final_score


def print_alignment(seq1, seq2, aligned_seq1, aligned_seq2, score):
    print("Original Sequences:")
    print(f"S1: {seq1}")
    print(f"S2: {seq2}")
    print()

    print("Aligned Sequences:")
    print(f"S1: {aligned_seq1}")

    match_string = ""
    for i in range(len(aligned_seq1)):
        if aligned_seq1[i] == aligned_seq2[i]:
            match_string += "|"
        elif aligned_seq1[i] == "-" or aligned_seq2[i] == "-":
            match_string += " "
        else:
            match_string += "x"
    print(f"    {match_string}")

    print(f"S2: {aligned_seq2}")
    print()
    print(f"Alignment Score: {score}")


def print_score_matrix(matrix, seq1, seq2):
    print("\nScoring Matrix:")
    print("      -", end="")
    for char in seq2:
        print(f"    {char}", end="")
    print()

    for i in range(len(matrix)):
        if i == 0:
            print(f"  -", end="")
        else:
            print(f"  {seq1[i-1]}", end="")

        for j in range(len(matrix[i])):
            print(f" {matrix[i][j]:4}", end="")
        print()


if __name__ == "__main__":
    S1 = "ACCGTGAAGCCAATAC"
    S2 = "AGCGTGCAGCCAATAC"

    print("=" * 80)
    print("NEEDLEMAN-WUNSCH GLOBAL ALIGNMENT ALGORITHM")
    print("=" * 80)
    print()

    aligned_s1, aligned_s2, score_matrix, final_score = needleman_wunsch(S1, S2)

    print_alignment(S1, S2, aligned_s1, aligned_s2, final_score)
    print_score_matrix(score_matrix, S1, S2)

    matches = sum(1 for i in range(len(aligned_s1)) if aligned_s1[i] == aligned_s2[i])
    mismatches = sum(1 for i in range(len(aligned_s1))
                    if aligned_s1[i] != aligned_s2[i] and aligned_s1[i] != "-" and aligned_s2[i] != "-")
    gaps = sum(1 for i in range(len(aligned_s1)) if aligned_s1[i] == "-" or aligned_s2[i] == "-")

    print("\n" + "=" * 80)
    print("ALIGNMENT STATISTICS")
    print("=" * 80)
    print(f"Alignment Length: {len(aligned_s1)}")
    print(f"Matches: {matches}")
    print(f"Mismatches: {mismatches}")
    print(f"Gaps: {gaps}")
    print(f"Identity: {matches / len(aligned_s1) * 100:.2f}%")
