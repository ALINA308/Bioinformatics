"""
Exercise 3: Three Scoring Equations for Sequence Similarity

1. Percent Identity (PID): Measures the percentage of exact matches
2. Similarity Score with Transition/Transversion Weighting (STW)
3. Normalized Alignment Score (NAS): Normalizes the raw alignment score
"""


def percent_identity_score(aligned_seq1, aligned_seq2):
    """Formula: PID = (M / L) * 100"""
    if len(aligned_seq1) != len(aligned_seq2):
        raise ValueError("Aligned sequences must have the same length")

    length = len(aligned_seq1)
    matches = sum(1 for i in range(length) if aligned_seq1[i] == aligned_seq2[i])

    pid = (matches / length) * 100 if length > 0 else 0

    return {
        'score_name': 'Percent Identity (PID)',
        'formula': 'PID = (M / L) * 100',
        'score': pid,
        'matches': matches,
        'alignment_length': length,
        'description': 'Percentage of exact matches in the alignment'
    }


def similarity_score_transition_transversion(aligned_seq1, aligned_seq2):
    """Formula: STW = (M * w_match + Tr * w_transition + Tv * w_transversion + G * w_gap) / L"""
    if len(aligned_seq1) != len(aligned_seq2):
        raise ValueError("Aligned sequences must have the same length")

    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}

    w_match = 1.0
    w_transition = 0.5
    w_transversion = 0.25
    w_gap = 0.0

    length = len(aligned_seq1)
    matches = 0
    transitions = 0
    transversions = 0
    gaps = 0

    for i in range(length):
        base1 = aligned_seq1[i]
        base2 = aligned_seq2[i]

        if base1 == base2:
            matches += 1
        elif base1 == '-' or base2 == '-':
            gaps += 1
        elif (base1 in purines and base2 in purines) or (base1 in pyrimidines and base2 in pyrimidines):
            transitions += 1
        else:
            transversions += 1

    numerator = (matches * w_match +
                transitions * w_transition +
                transversions * w_transversion +
                gaps * w_gap)

    stw = (numerator / length) * 100 if length > 0 else 0

    return {
        'score_name': 'Similarity Score with Transition/Transversion Weighting (STW)',
        'formula': 'STW = (M*w_m + Tr*w_tr + Tv*w_tv + G*w_g) / L * 100',
        'score': stw,
        'matches': matches,
        'transitions': transitions,
        'transversions': transversions,
        'gaps': gaps,
        'alignment_length': length,
        'weights': {
            'match': w_match,
            'transition': w_transition,
            'transversion': w_transversion,
            'gap': w_gap
        },
        'description': 'Weighted similarity score accounting for biological mutation patterns'
    }


def normalized_alignment_score(aligned_seq1, aligned_seq2, match_score=2, mismatch_penalty=-1, gap_penalty=-2):
    """Formula: NAS = (S - S_random) / (S_perfect - S_random) * 100"""
    if len(aligned_seq1) != len(aligned_seq2):
        raise ValueError("Aligned sequences must have the same length")

    length = len(aligned_seq1)

    raw_score = 0
    for i in range(length):
        base1 = aligned_seq1[i]
        base2 = aligned_seq2[i]

        if base1 == base2 and base1 != '-':
            raw_score += match_score
        elif base1 == '-' or base2 == '-':
            raw_score += gap_penalty
        else:
            raw_score += mismatch_penalty

    perfect_score = length * match_score

    expected_match_prob = 0.25
    expected_mismatch_prob = 0.75

    random_score = length * (match_score * expected_match_prob + mismatch_penalty * expected_mismatch_prob)

    if perfect_score != random_score:
        nas = ((raw_score - random_score) / (perfect_score - random_score)) * 100
    else:
        nas = 0

    nas = max(0, min(100, nas))

    return {
        'score_name': 'Normalized Alignment Score (NAS)',
        'formula': 'NAS = (S - S_random) / (S_perfect - S_random) * 100',
        'score': nas,
        'raw_score': raw_score,
        'perfect_score': perfect_score,
        'random_score': random_score,
        'alignment_length': length,
        'scoring_params': {
            'match_score': match_score,
            'mismatch_penalty': mismatch_penalty,
            'gap_penalty': gap_penalty
        },
        'description': 'Normalized score showing how much better than random the alignment is'
    }


def calculate_all_scores(aligned_seq1, aligned_seq2, match_score=2, mismatch_penalty=-1, gap_penalty=-2):
    score1 = percent_identity_score(aligned_seq1, aligned_seq2)
    score2 = similarity_score_transition_transversion(aligned_seq1, aligned_seq2)
    score3 = normalized_alignment_score(aligned_seq1, aligned_seq2, match_score, mismatch_penalty, gap_penalty)

    return {
        'pid': score1,
        'stw': score2,
        'nas': score3
    }


def print_score_comparison(scores):
    print("\n" + "=" * 100)
    print("SCORING EQUATIONS COMPARISON")
    print("=" * 100)

    for key, score in scores.items():
        print(f"\n{score['score_name']}")
        print("-" * 100)
        print(f"Formula: {score['formula']}")
        print(f"Score: {score['score']:.2f}")
        print(f"Description: {score['description']}")

        if key == 'pid':
            print(f"Details: {score['matches']} matches out of {score['alignment_length']} positions")
        elif key == 'stw':
            print(f"Details:")
            print(f"  - Matches: {score['matches']}")
            print(f"  - Transitions: {score['transitions']}")
            print(f"  - Transversions: {score['transversions']}")
            print(f"  - Gaps: {score['gaps']}")
            print(f"  - Weights: match={score['weights']['match']}, "
                  f"transition={score['weights']['transition']}, "
                  f"transversion={score['weights']['transversion']}, "
                  f"gap={score['weights']['gap']}")
        elif key == 'nas':
            print(f"Details:")
            print(f"  - Raw score: {score['raw_score']}")
            print(f"  - Perfect score: {score['perfect_score']}")
            print(f"  - Random score: {score['random_score']:.2f}")


def parse_fasta(filename):
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                sequence += line.upper()
    return sequence


def apply_scores_to_actual_alignment():
    print("\n" + "=" * 100)
    print("APPLYING SCORING EQUATIONS TO ACTUAL GENOME ALIGNMENT")
    print("=" * 100)

    from ex2_genome_alignment import smith_waterman, traceback_smith_waterman

    print("\nLoading genomes...")
    try:
        influenza_seq = parse_fasta("influenza_genome.fasta")
        covid_seq = parse_fasta("covid19_genome.fasta")

        print(f"Influenza: {len(influenza_seq)} bp")
        print(f"COVID-19: {len(covid_seq)} bp")

        chunk_size = 500
        chunk = influenza_seq[:chunk_size]

        print(f"\nAligning first {chunk_size} bp of Influenza against COVID-19 genome...")

        max_score, max_pos, score_matrix = smith_waterman(chunk, covid_seq)
        aligned_seq1, aligned_seq2, start_pos = traceback_smith_waterman(
            chunk, covid_seq, score_matrix, max_pos
        )

        print(f"Alignment found: length {len(aligned_seq1)} bp")

        scores = calculate_all_scores(aligned_seq1, aligned_seq2)

        print("\nAlignment preview (first 100 positions):")
        print(f"Influenza: {aligned_seq1[:100]}")
        match_str = ""
        for i in range(min(100, len(aligned_seq1))):
            if aligned_seq1[i] == aligned_seq2[i]:
                match_str += "|"
            elif aligned_seq1[i] == "-" or aligned_seq2[i] == "-":
                match_str += " "
            else:
                match_str += "x"
        print(f"           {match_str}")
        print(f"COVID-19:  {aligned_seq2[:100]}")

        print_score_comparison(scores)

        print("\n" + "=" * 100)
        print("SCORE INTERPRETATION")
        print("=" * 100)
        print(f"\nPID = {scores['pid']['score']:.2f}%")
        print("  Interpretation: Raw percentage of identical bases in the alignment")
        print(f"\nSTW = {scores['stw']['score']:.2f}%")
        print("  Interpretation: Weighted similarity considering biological mutation patterns")
        print("  (Transitions are weighted higher than transversions)")
        print(f"\nNAS = {scores['nas']['score']:.2f}%")
        print("  Interpretation: Normalized score showing how much better than random")
        print("  (0% = random, 100% = perfect match)")

    except FileNotFoundError as e:
        print(f"\nError: Could not find genome files. Please run ex2_genome_alignment.py first.")
        print(f"Details: {e}")


if __name__ == "__main__":
    print("=" * 100)
    print("THREE SCORING EQUATIONS FOR SEQUENCE SIMILARITY")
    print("=" * 100)

    print("\n" + "=" * 100)
    print("DEMONSTRATION WITH EXAMPLE SEQUENCES")
    print("=" * 100)

    print("\nExample 1: DNA sequences from Exercise 1")
    seq1 = "ACCGTGAAGCCAATAC"
    seq2 = "AGCGTGCAGCCAATAC"
    print(f"Seq1: {seq1}")
    print(f"Seq2: {seq2}")

    scores = calculate_all_scores(seq1, seq2)
    print_score_comparison(scores)

    print("\n" + "=" * 100)
    print("\nExample 2: Sequences with gaps")
    seq3 = "ACGT-GCAT-CGA"
    seq4 = "A-GTTG-ATACGA"
    print(f"Seq1: {seq3}")
    print(f"Seq2: {seq4}")

    scores2 = calculate_all_scores(seq3, seq4)
    print_score_comparison(scores2)

    print("\n" + "=" * 100)
    print("\nExample 3: Perfect match")
    seq5 = "AAATTTGGGCCC"
    seq6 = "AAATTTGGGCCC"
    print(f"Seq1: {seq5}")
    print(f"Seq2: {seq6}")

    scores3 = calculate_all_scores(seq5, seq6)
    print_score_comparison(scores3)

    apply_scores_to_actual_alignment()

    print("\n" + "=" * 100)
    print("SUMMARY OF SCORING EQUATIONS")
    print("=" * 100)
    print("\n1. Percent Identity (PID):")
    print("   - Simple and intuitive")
    print("   - Counts exact matches only")
    print("   - Range: 0-100%")
    print("   - Best for: Quick similarity assessment")

    print("\n2. Similarity Score with Transition/Transversion Weighting (STW):")
    print("   - Biologically informed")
    print("   - Accounts for mutation patterns")
    print("   - Transitions weighted higher than transversions")
    print("   - Best for: Evolutionary distance estimation")

    print("\n3. Normalized Alignment Score (NAS):")
    print("   - Statistically rigorous")
    print("   - Accounts for expected random similarity")
    print("   - Normalized to 0-100% scale")
    print("   - Best for: Statistical significance assessment")

    print("\n" + "=" * 100)
