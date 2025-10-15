
import math

def calculate_tm_simple(dna_sequence):
    
    dna_sequence = dna_sequence.upper()

    g_count = dna_sequence.count('G')
    c_count = dna_sequence.count('C')
    a_count = dna_sequence.count('A')
    t_count = dna_sequence.count('T')

    tm = 4 * (g_count + c_count) + 2 * (a_count + t_count)

    return tm


def calculate_tm_salt_adjusted(dna_sequence, na_concentration=50):
   
    dna_sequence = dna_sequence.upper()
    length = len(dna_sequence)

    if length == 0:
        return 0

    g_count = dna_sequence.count('G')
    c_count = dna_sequence.count('C')

    gc_content = ((g_count + c_count) / length) * 100

   
    na_molar = na_concentration / 1000

    tm = 81.5 + 16.6 * math.log10(na_molar) + 0.41 * gc_content - (600 / length)

    return tm


def validate_dna_sequence(dna_sequence):
    
    valid_nucleotides = set('ATGCatgc')
    return all(nucleotide in valid_nucleotides for nucleotide in dna_sequence)


def main():
    print("=" * 60)
    print("DNA Melting Temperature (Tm) Calculator")
    print("=" * 60)
    print()

   
    dna_sequence = input("Enter DNA sequence: ").strip()

   
    if not dna_sequence:
        print("Error: Empty sequence provided.")
        return

    if not validate_dna_sequence(dna_sequence):
        print("Error: Invalid DNA sequence. Only A, T, G, C nucleotides are allowed.")
        return

    print()
    print("-" * 60)
    print(f"DNA Sequence: {dna_sequence.upper()}")
    print(f"Length: {len(dna_sequence)} bp")
    print("-" * 60)
    print()

   
    tm_simple = calculate_tm_simple(dna_sequence)
    print("Method 1: Simple Formula")
    print("Formula: Tm = 4(G + C) + 2(A + T)")
    print(f"Result: {tm_simple:.2f} C")
    print()

   
    print("Method 2: Salt-Adjusted Formula")
    print("Formula: Tm = 81.5 + 16.6(log10([Na+])) + 0.41*(%GC) - 600/length")

    
    use_default = input("Use default Na+ concentration (50 mM)? (y/n): ").strip().lower()

    if use_default == 'n':
        try:
            na_conc = float(input("Enter Na+ concentration in mM: "))
            if na_conc <= 0:
                print("Warning: Invalid concentration. Using default 50 mM.")
                na_conc = 50
        except ValueError:
            print("Warning: Invalid input. Using default 50 mM.")
            na_conc = 50
    else:
        na_conc = 50

    tm_salt = calculate_tm_salt_adjusted(dna_sequence, na_conc)
    print(f"Na+ concentration: {na_conc} mM")
    print(f"Result: {tm_salt:.2f} C")
    print()
    print("=" * 60)


if __name__ == "__main__":
    main()
