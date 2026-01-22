# -*- coding: utf-8 -*-
import numpy as np
from collections import defaultdict
import re

# Sample poetry texts (excerpts for demonstration)
# Mihai Eminescu - classic, romantic style with nature imagery
eminescu_text = """
Somnoroase pasarele pe la cuiburi se aduna.
Se ascund in ramurele si dorm in liniste.
Noapte buna si pace pentru toate.

Doar izvoarele suspina cand codru negru tace.
Dorm si florile in gradina sub luna care straluceste.
Dormi in pace si somn dulce sa ai.

Trece lebada pe ape intre trestii sa se culce.
Fie ingerii aproape in noapte cu somnul dulce.
Luna straluceste peste dealuri si peste vai.

Peste noptii se ridica mandra luna care straluceste.
Totul este vis si armonie in noapte.
Noapte buna si somn dulce pentru tine.

Luna este mare si frumoasa in cer.
Straluceste peste codru si peste vale.
Noapte vine cu pace si somn pentru toate.

Florile dorm in gradina sub luna.
Pasarele se ascund in ramuri si dorm.
Totul este liniste si pace in noapte.

Apa curge incet si lebada doarme.
Luna straluceste si totul este frumos.
Noapte buna si somn dulce sa fie.
"""

# Nichita Stanescu - modern, abstract style with existential themes
stanescu_text = """
Sunt un om foarte singur asa cum e singur primul om.
Sunt singur si nu stiu ce simt in mine.
Poate ca sunt un gol foarte mare.

Nu stiu sa vorbesc despre asta.
Nu stiu ce este in mine.
Sunt singur asa cum va fi singur ultimul om.

Poate ca nu e nimeni in mine si sunt gol.
Poate ca sunt un ochi care priveste fara sa vada.
Nu stiu ce simt si nu stiu ce este.

Sunt un om foarte ciudat care nu stie ce este el.
Sunt un om care cauta in locuri unde nu este nimic.
Sunt singur si totusi nu simt nici macar tristetea.

Totul este gol in mine si nu stiu.
Nu stiu ce caut si nu stiu ce simt.
Sunt un om singur foarte singur.

Poate ca totul este un gol foarte mare.
Nu stiu sa vorbesc si nu stiu sa simt.
Sunt singur asa cum e singur omul.

Cauta in mine si nu este nimic.
Totul este gol si foarte singur.
Nu stiu ce este si nu stiu ce simt.
"""

# Simulated mixed text (created by combining styles)
# This simulates what "Mihai" (the accused) might have written
# Part 1: Eminescu-like (nature, romantic)
# Part 2: Stanescu-like (existential, abstract)
# Part 3: Mixed style
mihai_text = """
Luna straluceste peste dealuri si peste vai.
Noapte vine cu pace si somn pentru toate.
Florile dorm in gradina sub luna care straluceste.
Totul este liniste si pace in noapte.

Sunt un om foarte singur si nu stiu ce simt.
Poate ca sunt un gol foarte mare in mine.
Nu stiu ce este si nu stiu sa vorbesc.
Sunt singur asa cum e singur omul.

Luna este mare si totul este gol.
Poate ca noapte vine si sunt singur.
Nu stiu ce simt in noapte sub luna.
Totul este liniste si totul este gol.
"""


def preprocess_text(text):
    """Convert text to lowercase and split into words."""
    # Remove newlines and extra spaces
    text = text.lower()
    # Keep only letters and spaces (simple approach without diacritics)
    text = re.sub(r'[^a-z\s]', ' ', text)
    # Split into words and remove empty strings
    words = [w for w in text.split() if w]
    return words


def build_transition_matrix(words, epsilon=1e-5):
    """
    Build word transition probability matrix.

    Args:
        words: List of words
        epsilon: Pseudocount to avoid log(0)

    Returns:
        vocab: Vocabulary (unique words)
        prob_matrix: Dictionary of transition probabilities
    """
    # Build vocabulary
    vocab = set(words)

    # Count transitions (word -> next_word)
    transition_counts = defaultdict(lambda: defaultdict(float))

    # Initialize with epsilon
    for w1 in vocab:
        for w2 in vocab:
            transition_counts[w1][w2] = epsilon

    # Count actual transitions
    for i in range(len(words) - 1):
        transition_counts[words[i]][words[i+1]] += 1.0

    # Convert to probabilities
    prob_matrix = {}
    for w1 in vocab:
        total = sum(transition_counts[w1].values())
        prob_matrix[w1] = {}
        for w2 in vocab:
            prob_matrix[w1][w2] = transition_counts[w1][w2] / total

    return vocab, prob_matrix


def score_window(words, ll_matrix, vocab_eminescu, vocab_stanescu):
    """
    Calculate log-likelihood score for a window of words.
    Positive = Eminescu-like, Negative = Stanescu-like
    """
    score = 0.0
    valid_transitions = 0

    for i in range(len(words) - 1):
        w1, w2 = words[i], words[i+1]

        # Check if both words exist in both vocabularies
        if w1 in vocab_eminescu and w2 in vocab_eminescu and \
           w1 in vocab_stanescu and w2 in vocab_stanescu:
            if w1 in ll_matrix and w2 in ll_matrix[w1]:
                score += ll_matrix[w1][w2]
                valid_transitions += 1

    # Return average score if we have valid transitions
    return score / valid_transitions if valid_transitions > 0 else 0.0


def main():
    print("=" * 80)
    print("PLAGIARISM DETECTION USING MARKOV MODELS")
    print("Courtroom Analysis: Identifying Author Styles")
    print("=" * 80)
    print()

    # Preprocess texts
    print("Step 1: Preprocessing texts...")
    eminescu_words = preprocess_text(eminescu_text)
    stanescu_words = preprocess_text(stanescu_text)
    mihai_words = preprocess_text(mihai_text)

    print(f"  Eminescu: {len(eminescu_words)} words")
    print(f"  Stanescu: {len(stanescu_words)} words")
    print(f"  Accused text: {len(mihai_words)} words")
    print()

    # Build transition matrices
    print("Step 2: Building transition probability matrices...")
    vocab_eminescu, prob_eminescu = build_transition_matrix(eminescu_words)
    vocab_stanescu, prob_stanescu = build_transition_matrix(stanescu_words)
    print(f"  Eminescu vocabulary: {len(vocab_eminescu)} unique words")
    print(f"  Stanescu vocabulary: {len(vocab_stanescu)} unique words")
    print()

    # Build log-likelihood matrix
    print("Step 3: Creating log-likelihood matrix...")
    # Common vocabulary
    common_vocab = vocab_eminescu & vocab_stanescu
    print(f"  Common vocabulary: {len(common_vocab)} words")

    ll_matrix = {}
    for w1 in common_vocab:
        ll_matrix[w1] = {}
        for w2 in common_vocab:
            # Log-likelihood ratio: log2(P_eminescu / P_stanescu)
            # Positive values = more like Eminescu
            # Negative values = more like Stanescu
            ll_matrix[w1][w2] = np.log2(
                prob_eminescu[w1][w2] / prob_stanescu[w1][w2]
            )
    print()

    # Analyze the accused text using sliding window
    print("Step 4: Analyzing accused text with sliding window...")
    print("=" * 80)

    window_size = 10  # Number of words per window
    threshold = 0.5   # Threshold for classification

    results = []

    for i in range(0, len(mihai_words) - window_size + 1, 3):  # Step by 3 words
        window = mihai_words[i:i+window_size]
        score = score_window(window, ll_matrix, vocab_eminescu, vocab_stanescu)

        # Classify based on score
        if score > threshold:
            author = "EMINESCU"
        elif score < -threshold:
            author = "STANESCU"
        else:
            author = "UNCERTAIN/MIXED"

        results.append({
            'position': i,
            'window': ' '.join(window[:5]) + '...',
            'score': score,
            'author': author
        })

    # Display results
    print(f"{'Position':<10} {'Window Preview':<30} {'Score':<10} {'Attribution':<20}")
    print("-" * 80)

    for result in results:
        print(f"{result['position']:<10} {result['window']:<30} "
              f"{result['score']:>8.3f}  {result['author']:<20}")

    print()
    print("=" * 80)
    print("INTERPRETATION:")
    print("  Score > 0.5  -> Text resembles MIHAI EMINESCU's style (classic, romantic)")
    print("  Score < -0.5 -> Text resembles NICHITA STANESCU's style (modern, abstract)")
    print("  -0.5 <= Score <= 0.5 -> Mixed or uncertain attribution")
    print("=" * 80)
    print()

    # Summary statistics
    eminescu_count = sum(1 for r in results if r['author'] == 'EMINESCU')
    stanescu_count = sum(1 for r in results if r['author'] == 'STANESCU')
    uncertain_count = sum(1 for r in results if r['author'] == 'UNCERTAIN/MIXED')

    print("VERDICT SUMMARY:")
    print(f"  Windows attributed to Eminescu: {eminescu_count}")
    print(f"  Windows attributed to Stanescu: {stanescu_count}")
    print(f"  Uncertain/Mixed windows: {uncertain_count}")
    print()

    if eminescu_count > 0 and stanescu_count > 0:
        print("CONCLUSION: The accused text contains passages from BOTH authors.")
        print("Evidence of PLAGIARISM detected!")
    elif eminescu_count > stanescu_count:
        print("CONCLUSION: The text is primarily in Eminescu's style.")
    elif stanescu_count > eminescu_count:
        print("CONCLUSION: The text is primarily in Stanescu's style.")
    else:
        print("CONCLUSION: Attribution is uncertain. More data needed.")
    print("=" * 80)


if __name__ == "__main__":
    main()
