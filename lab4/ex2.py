"""
 Download COVID-19 and Influenza genomes from NCBI,
compare codon frequencies, create charts, and analyze amino acids.
"""
from urllib.request import urlopen
from collections import Counter
import matplotlib.pyplot as plt
from ex1 import _normalize_seq, translate_coding_region, STOP_CODONS, CODON_TABLE

def download_fasta(url, filename):
    
    try:
        with urlopen(url) as response:
            data = response.read().decode('utf-8')
        with open(filename, 'w') as f:
            f.write(data)
        print(f"Downloaded {filename}")
    except Exception as e:
        print(f"Error downloading {filename}: {e}")

def parse_fasta(filename):
   
    with open(filename) as f:
        lines = f.readlines()
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return seq

def count_codons(seq):
   
    rna = _normalize_seq(seq)
    codon_count = Counter()
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        codon_count[codon] += 1
    return codon_count

def plot_top_codons(codon_count, title, filename):

    top = sorted(codon_count.items(), key=lambda x: x[1], reverse=True)[:10]
    codons, counts = zip(*top)
    plt.figure(figsize=(10, 6))
    plt.bar(codons, counts)
    plt.title(title)
    plt.xticks(rotation=45)
    plt.ylabel('Frequency')
    plt.savefig(filename)
    plt.close()

def plot_comparison(covid_codons, flu_codons):
    top_covid = dict(sorted(covid_codons.items(), key=lambda x: x[1], reverse=True)[:10])
    top_flu = dict(sorted(flu_codons.items(), key=lambda x: x[1], reverse=True)[:10])
    all_codons = set(top_covid.keys()) | set(top_flu.keys())
    covid_vals = [top_covid.get(c, 0) for c in all_codons]
    flu_vals = [top_flu.get(c, 0) for c in all_codons]
    x = range(len(all_codons))
    plt.figure(figsize=(12, 6))
    plt.bar(x, covid_vals, width=0.4, label='COVID-19', align='center')
    plt.bar([i + 0.4 for i in x], flu_vals, width=0.4, label='Influenza', align='center')
    plt.xticks([i + 0.2 for i in x], list(all_codons), rotation=45)
    plt.legend()
    plt.title('Comparison of Top Codons')
    plt.ylabel('Frequency')
    plt.savefig('comparison_codons.png')
    plt.close()

def print_top_codons(codon_count, name):
    top = sorted(codon_count.items(), key=lambda x: x[1], reverse=True)[:10]
    codons = [codon for codon, _ in top]
    print(f"Top 10 codons for {name}: {codons}")

def print_top_aa(aa_count, name):
    top = sorted(aa_count.items(), key=lambda x: x[1], reverse=True)[:3]
    print(f"Top 3 Amino Acids in {name}: {top}")

def main():
    
    covid_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.1&rettype=fasta&retmode=text"
    flu_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_002023.1&rettype=fasta&retmode=text"
    
    download_fasta(covid_url, 'covid19.fasta')
    download_fasta(flu_url, 'influenza.fasta')

    covid_seq = parse_fasta('covid19.fasta')
    flu_seq = parse_fasta('influenza.fasta')

    covid_codons = count_codons(covid_seq)
    flu_codons = count_codons(flu_seq)

    # Translate and count amino acids
    covid_protein, _ = translate_coding_region(covid_seq)
    flu_protein, _ = translate_coding_region(flu_seq)
    covid_aa = Counter(covid_protein)
    flu_aa = Counter(flu_protein)

    # Count amino acids from codons
    covid_aa_from_codons = Counter()
    for codon, count in covid_codons.items():
        aa = CODON_TABLE.get(codon, 'X')
        covid_aa_from_codons[aa] += count
    flu_aa_from_codons = Counter()
    for codon, count in flu_codons.items():
        aa = CODON_TABLE.get(codon, 'X')
        flu_aa_from_codons[aa] += count

   
    print_top_codons(covid_codons, 'COVID-19')
    print_top_codons(flu_codons, 'Influenza')
    plot_top_codons(covid_codons, 'Top 10 Codons COVID-19', 'covid_codons.png')
    plot_top_codons(flu_codons, 'Top 10 Codons Influenza', 'flu_codons.png')
    plot_comparison(covid_codons, flu_codons)
    print_top_aa(covid_aa_from_codons, 'SARS-CoV')
    print_top_aa(flu_aa_from_codons, 'Influenza A')

    
if __name__ == '__main__':
    main()
