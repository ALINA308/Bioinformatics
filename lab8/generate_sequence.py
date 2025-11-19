import random
import json

def generate_random_dna(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[base] for base in reversed(seq))

def create_transposon(length, tir_length=10, tsd_length=8):
    tsd = generate_random_dna(tsd_length)
    tir = generate_random_dna(tir_length)
    tir_rc = reverse_complement(tir)
    internal_length = length - 2 * tir_length - 2 * tsd_length
    if internal_length < 10:
        internal_length = 10
    internal_seq = generate_random_dna(internal_length)
    transposon = tsd + tir + internal_seq + tir_rc + tsd

    return {
        'sequence': transposon,
        'tsd': tsd,
        'tir': tir,
        'tir_rc': tir_rc,
        'internal': internal_seq,
        'total_length': len(transposon)
    }

def insert_transposons_into_sequence(sequence_length, num_transposons=3):
    transposon_lengths = []
    min_transposon_size = 40
    max_transposon_size = 80

    for _ in range(num_transposons):
        transposon_lengths.append(random.randint(min_transposon_size, max_transposon_size))

    total_transposon_length = sum(transposon_lengths)
    remaining_length = sequence_length - total_transposon_length

    if remaining_length < num_transposons + 1:
        total_transposon_length = int(sequence_length * 0.7)
        avg_length = total_transposon_length // num_transposons
        transposon_lengths = [avg_length] * num_transposons
        remaining_length = sequence_length - total_transposon_length

    spacer_lengths = []
    for _ in range(num_transposons + 1):
        spacer_lengths.append(max(5, remaining_length // (num_transposons + 1)))

    sequence = ""
    transposons_info = []
    current_pos = 0

    for i in range(num_transposons):
        spacer = generate_random_dna(spacer_lengths[i])
        sequence += spacer
        current_pos += len(spacer)

        transposon = create_transposon(transposon_lengths[i])
        start_pos = current_pos
        sequence += transposon['sequence']
        end_pos = current_pos + transposon['total_length'] - 1

        transposons_info.append({
            'id': i + 1,
            'start': start_pos,
            'end': end_pos,
            'length': transposon['total_length'],
            'tsd': transposon['tsd'],
            'tir': transposon['tir'],
            'tir_rc': transposon['tir_rc'],
            'sequence': transposon['sequence']
        })

        current_pos += transposon['total_length']

    final_spacer = generate_random_dna(spacer_lengths[-1])
    sequence += final_spacer

    return {
        'sequence': sequence,
        'length': len(sequence),
        'transposons': transposons_info,
        'num_transposons': num_transposons
    }

def main():
    random.seed(42)
    num_transposons = random.randint(3, 4)
    sequence_length = random.randint(200, 400)

    print(f"Generating DNA sequence of length {sequence_length}bp with {num_transposons} transposable elements...")

    result = insert_transposons_into_sequence(sequence_length, num_transposons)

    with open('artificial_dna.fasta', 'w') as f:
        f.write(f">Artificial_DNA_Sequence length={result['length']}\n")
        seq = result['sequence']
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + '\n')

    with open('transposons_ground_truth.json', 'w') as f:
        json.dump(result, f, indent=2)

    print(f"\n{'='*70}")
    print(f"GENERATED DNA SEQUENCE")
    print(f"{'='*70}")
    print(f"Total length: {result['length']} bp")
    print(f"Number of transposons: {result['num_transposons']}")
    print(f"\nTransposon positions (ground truth):")
    print(f"{'ID':<5} {'Start':<8} {'End':<8} {'Length':<8} {'TSD':<12} {'TIR':<15}")
    print(f"{'-'*70}")

    for te in result['transposons']:
        print(f"{te['id']:<5} {te['start']:<8} {te['end']:<8} {te['length']:<8} {te['tsd']:<12} {te['tir']:<15}")

    print(f"\nFiles created:")
    print(f"  - artificial_dna.fasta (DNA sequence)")
    print(f"  - transposons_ground_truth.json (true positions)")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
