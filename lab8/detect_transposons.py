import json
import re
from collections import defaultdict

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[base] for base in reversed(seq))

def find_inverted_repeats(sequence, min_length=8, max_length=15, max_gap=20, max_distance=100):
    inverted_repeats = []

    for i in range(len(sequence) - min_length):
        for length in range(min_length, min(max_length + 1, len(sequence) - i + 1)):
            left_seq = sequence[i:i+length]

            for j in range(i + length + max_gap, min(i + max_distance, len(sequence) - length + 1)):
                right_seq = sequence[j:j+length]

                if left_seq == reverse_complement(right_seq):
                    inverted_repeats.append({
                        'left_start': i,
                        'left_end': i + length - 1,
                        'right_start': j,
                        'right_end': j + length - 1,
                        'repeat_seq': left_seq,
                        'length': length,
                        'internal_size': j - (i + length)
                    })

    filtered_repeats = []
    inverted_repeats.sort(key=lambda x: -x['length'])

    used_positions = set()
    for ir in inverted_repeats:
        positions = set(range(ir['left_start'], ir['left_end'] + 1))
        positions.update(range(ir['right_start'], ir['right_end'] + 1))

        if not positions.intersection(used_positions):
            filtered_repeats.append(ir)
            used_positions.update(positions)

    return filtered_repeats

def find_direct_repeats(sequence, min_length=5, max_length=12, max_distance=100):
    direct_repeats = []

    for i in range(len(sequence) - min_length):
        for length in range(min_length, min(max_length + 1, len(sequence) - i + 1)):
            repeat = sequence[i:i+length]

            for j in range(i + length, min(i + max_distance, len(sequence) - length + 1)):
                if sequence[j:j+length] == repeat:
                    distance = j - i
                    direct_repeats.append({
                        'first_start': i,
                        'first_end': i + length - 1,
                        'second_start': j,
                        'second_end': j + length - 1,
                        'repeat_seq': repeat,
                        'length': length,
                        'distance': distance
                    })

    return direct_repeats

def detect_transposons(sequence, min_te_length=30):
    print("Searching for inverted repeats (TIRs)...")
    inverted_repeats = find_inverted_repeats(sequence)
    print(f"Found {len(inverted_repeats)} inverted repeat pairs")

    print("\nSearching for direct repeats (TSDs)...")
    direct_repeats = find_direct_repeats(sequence)
    print(f"Found {len(direct_repeats)} direct repeat pairs")

    transposons = []

    for ir in inverted_repeats:
        tsd_start = ir['left_start']
        tsd_end = ir['right_end']

        matching_tsds = []
        for dr in direct_repeats:
            if (dr['first_end'] < ir['left_start'] and
                dr['second_start'] > ir['right_end'] and
                dr['second_start'] - dr['first_end'] > min_te_length):

                matching_tsds.append(dr)

        if matching_tsds:
            best_tsd = max(matching_tsds, key=lambda x: x['length'])

            transposons.append({
                'start': best_tsd['first_start'],
                'end': best_tsd['second_end'],
                'length': best_tsd['second_end'] - best_tsd['first_start'] + 1,
                'tir_left': ir['repeat_seq'],
                'tir_right': reverse_complement(ir['repeat_seq']),
                'tir_length': ir['length'],
                'tir_left_pos': (ir['left_start'], ir['left_end']),
                'tir_right_pos': (ir['right_start'], ir['right_end']),
                'tsd': best_tsd['repeat_seq'],
                'tsd_length': best_tsd['length'],
                'tsd_left_pos': (best_tsd['first_start'], best_tsd['first_end']),
                'tsd_right_pos': (best_tsd['second_start'], best_tsd['second_end']),
                'confidence': 'high'
            })
        else:
            transposons.append({
                'start': ir['left_start'],
                'end': ir['right_end'],
                'length': ir['right_end'] - ir['left_start'] + 1,
                'tir_left': ir['repeat_seq'],
                'tir_right': reverse_complement(ir['repeat_seq']),
                'tir_length': ir['length'],
                'tir_left_pos': (ir['left_start'], ir['left_end']),
                'tir_right_pos': (ir['right_start'], ir['right_end']),
                'tsd': None,
                'tsd_length': 0,
                'confidence': 'medium'
            })

    transposons.sort(key=lambda x: x['start'])
    filtered_transposons = []
    for te in transposons:
        overlap = False
        for existing_te in filtered_transposons:
            overlap_start = max(te['start'], existing_te['start'])
            overlap_end = min(te['end'], existing_te['end'])
            if overlap_end >= overlap_start:
                overlap = True
                break

        if not overlap:
            filtered_transposons.append(te)

    return filtered_transposons

def read_fasta(filename):
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

def main():
    print("="*70)
    print("TRANSPOSABLE ELEMENT DETECTION")
    print("="*70)

    print("\nReading DNA sequence from artificial_dna.fasta...")
    try:
        sequence = read_fasta('artificial_dna.fasta')
        print(f"Sequence length: {len(sequence)} bp")
    except FileNotFoundError:
        print("Error: artificial_dna.fasta not found. Run generate_sequence.py first.")
        return

    print(f"\n{'-'*70}")

    detected_transposons = detect_transposons(sequence)

    print(f"\n{'-'*70}")
    print(f"DETECTION RESULTS")
    print(f"{'-'*70}")
    print(f"Detected {len(detected_transposons)} transposable element(s)\n")

    if detected_transposons:
        print(f"{'ID':<5} {'Start':<8} {'End':<8} {'Length':<8} {'TIR Len':<10} {'TSD':<12} {'Confidence':<12}")
        print(f"{'-'*70}")

        for idx, te in enumerate(detected_transposons, 1):
            tsd_display = te['tsd'] if te['tsd'] else 'N/A'
            print(f"{idx:<5} {te['start']:<8} {te['end']:<8} {te['length']:<8} "
                  f"{te['tir_length']:<10} {tsd_display:<12} {te['confidence']:<12}")

            print(f"      TIR: {te['tir_left']} ... {te['tir_right']}")
            if te['tsd']:
                print(f"      Structure: {te['tsd']}|{te['tir_left']}...{te['tir_right']}|{te['tsd']}")
            print()

    results = {
        'sequence_length': len(sequence),
        'num_detected': len(detected_transposons),
        'transposons': detected_transposons
    }

    with open('detection_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print(f"Results saved to detection_results.json")

    try:
        with open('transposons_ground_truth.json', 'r') as f:
            ground_truth = json.load(f)

        print(f"\n{'-'*70}")
        print(f"COMPARISON WITH GROUND TRUTH")
        print(f"{'-'*70}")
        print(f"Expected: {ground_truth['num_transposons']} transposons")
        print(f"Detected: {len(detected_transposons)} transposons")

        print(f"\nGround truth positions:")
        for te in ground_truth['transposons']:
            print(f"  TE {te['id']}: {te['start']}-{te['end']} (TSD: {te['tsd']}, TIR: {te['tir']})")

        true_positives = 0
        for detected in detected_transposons:
            for true_te in ground_truth['transposons']:
                if (detected['start'] <= true_te['end'] and detected['end'] >= true_te['start']):
                    start_diff = abs(detected['start'] - true_te['start'])
                    end_diff = abs(detected['end'] - true_te['end'])
                    if start_diff <= 10 and end_diff <= 10:
                        true_positives += 1
                        break

        precision = true_positives / len(detected_transposons) if detected_transposons else 0
        recall = true_positives / ground_truth['num_transposons']

        print(f"\nAccuracy metrics:")
        print(f"  True Positives: {true_positives}")
        print(f"  Precision: {precision:.2%}")
        print(f"  Recall: {recall:.2%}")

    except FileNotFoundError:
        print("\nGround truth file not found - skipping validation")

    print(f"\n{'='*70}\n")

if __name__ == "__main__":
    main()
