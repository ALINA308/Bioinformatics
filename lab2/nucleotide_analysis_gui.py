"""
Dinucleotide and Trinucleotide Analysis - generates all possible dinucleotide and trinucleotide combinations
and calculates their percentage occurrence in a given DNA sequence.
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox

def generate_dinucleotides():
    
    nucleotides = ['A', 'C', 'G', 'T']
    dinucleotides = []

    for n1 in nucleotides:
        for n2 in nucleotides:
            dinucleotides.append(n1 + n2)

    return dinucleotides

def generate_trinucleotides():
    
    nucleotides = ['A', 'C', 'G', 'T']
    trinucleotides = []

    for n1 in nucleotides:
        for n2 in nucleotides:
            for n3 in nucleotides:
                trinucleotides.append(n1 + n2 + n3)

    return trinucleotides

def count_occurrences(sequence, pattern):
    
    count = 0
    for i in range(len(sequence) - len(pattern) + 1):
        if sequence[i:i+len(pattern)] == pattern:
            count += 1
    return count

def calculate_percentage(sequence, pattern):
   
    pattern_length = len(pattern)
    total_possible = len(sequence) - pattern_length + 1

    if total_possible <= 0:
        return 0.0

    count = count_occurrences(sequence, pattern)
    percentage = (count / total_possible) * 100

    return percentage

class NucleotideAnalyzerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Dinucleotide and Trinucleotide Analysis")
        self.root.geometry("900x700")
        self.root.resizable(True, True)

        
        self.default_sequence = "ATTGTCCCAATCTGTTG"

        self.create_widgets()

    def create_widgets(self):
       
        title_frame = tk.Frame(self.root, bg="#2c3e50", pady=10)
        title_frame.pack(fill=tk.X)

        title_label = tk.Label(title_frame, text="NUCLEOTIDE ANALYSIS TOOL",
                              font=("Arial", 18, "bold"), bg="#2c3e50", fg="white")
        title_label.pack()

       
        input_frame = tk.LabelFrame(self.root, text="DNA Sequence Input",
                                   font=("Arial", 11, "bold"), padx=10, pady=10)
        input_frame.pack(fill=tk.X, padx=10, pady=10)

        tk.Label(input_frame, text="Enter DNA Sequence:", font=("Arial", 10)).pack(anchor=tk.W)

        self.sequence_entry = tk.Entry(input_frame, font=("Courier", 10), width=80)
        self.sequence_entry.pack(fill=tk.X, pady=5)
        self.sequence_entry.insert(0, self.default_sequence)

        
        button_frame = tk.Frame(input_frame)
        button_frame.pack(pady=5)

        self.analyze_btn = tk.Button(button_frame, text="Analyze Sequence",
                                     command=self.analyze_sequence,
                                     bg="#27ae60", fg="white", font=("Arial", 10, "bold"),
                                     padx=20, pady=5, cursor="hand2")
        self.analyze_btn.pack(side=tk.LEFT, padx=5)

        self.clear_btn = tk.Button(button_frame, text="Clear Results",
                                   command=self.clear_results,
                                   bg="#e74c3c", fg="white", font=("Arial", 10, "bold"),
                                   padx=20, pady=5, cursor="hand2")
        self.clear_btn.pack(side=tk.LEFT, padx=5)

       
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        
        self.dinuc_frame = tk.Frame(self.notebook)
        self.notebook.add(self.dinuc_frame, text="Dinucleotides (16 combinations)")

        self.dinuc_text = scrolledtext.ScrolledText(self.dinuc_frame,
                                                    font=("Courier", 10),
                                                    wrap=tk.WORD)
        self.dinuc_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        
        self.trinuc_frame = tk.Frame(self.notebook)
        self.notebook.add(self.trinuc_frame, text="Trinucleotides (64 combinations)")

        self.trinuc_text = scrolledtext.ScrolledText(self.trinuc_frame,
                                                     font=("Courier", 10),
                                                     wrap=tk.WORD)
        self.trinuc_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        
        self.summary_frame = tk.Frame(self.notebook)
        self.notebook.add(self.summary_frame, text="Summary")

        self.summary_text = scrolledtext.ScrolledText(self.summary_frame,
                                                      font=("Courier", 10),
                                                      wrap=tk.WORD)
        self.summary_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        
        self.status_bar = tk.Label(self.root, text="Ready to analyze",
                                  bd=1, relief=tk.SUNKEN, anchor=tk.W,
                                  font=("Arial", 9))
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)

    def clear_results(self):
       
        self.dinuc_text.delete(1.0, tk.END)
        self.trinuc_text.delete(1.0, tk.END)
        self.summary_text.delete(1.0, tk.END)
        self.status_bar.config(text="Results cleared")

    def analyze_sequence(self):
        
        sequence = self.sequence_entry.get().strip().upper()

        if not sequence:
            messagebox.showwarning("Input Error", "Please enter a DNA sequence!")
            return

       
        valid_nucleotides = set('ACGT')
        if not all(c in valid_nucleotides for c in sequence):
            messagebox.showerror("Invalid Sequence",
                               "Sequence must contain only A, C, G, T nucleotides!")
            return

        self.status_bar.config(text="Analyzing sequence...")
        self.root.update()

        
        self.clear_results()

        
        self.analyze_dinucleotides(sequence)

       
        self.analyze_trinucleotides(sequence)

        
        self.generate_summary(sequence)

        self.status_bar.config(text=f"Analysis complete! Sequence length: {len(sequence)} nucleotides")

    def analyze_dinucleotides(self, sequence):
       
        dinucleotides = generate_dinucleotides()

        output = "=" * 70 + "\n"
        output += "DINUCLEOTIDE ANALYSIS\n"
        output += "=" * 70 + "\n\n"
        output += f"Sequence: {sequence}\n"
        output += f"Sequence Length: {len(sequence)} nucleotides\n\n"
        output += f"Total possible dinucleotides: {len(dinucleotides)}\n"
        output += f"Total possible dinucleotide positions: {len(sequence) - 1}\n\n"

        # Calculate percentages
        dinuc_results = []
        for dinuc in dinucleotides:
            count = count_occurrences(sequence, dinuc)
            percentage = calculate_percentage(sequence, dinuc)
            dinuc_results.append((dinuc, count, percentage))

        # Sort by percentage (descending)
        dinuc_results.sort(key=lambda x: x[2], reverse=True)

        output += f"{'Dinucleotide':<15}{'Count':<10}{'Percentage':<15}\n"
        output += "-" * 40 + "\n"

        for dinuc, count, percentage in dinuc_results:
            output += f"{dinuc:<15}{count:<10}{percentage:>6.2f}%\n"

        self.dinuc_text.insert(1.0, output)
        self.dinuc_results = dinuc_results

    def analyze_trinucleotides(self, sequence):
       
        trinucleotides = generate_trinucleotides()

        output = "=" * 70 + "\n"
        output += "TRINUCLEOTIDE ANALYSIS\n"
        output += "=" * 70 + "\n\n"
        output += f"Sequence: {sequence}\n"
        output += f"Sequence Length: {len(sequence)} nucleotides\n\n"
        output += f"Total possible trinucleotides: {len(trinucleotides)}\n"
        output += f"Total possible trinucleotide positions: {len(sequence) - 2}\n\n"

        # Calculate percentages
        trinuc_results = []
        for trinuc in trinucleotides:
            count = count_occurrences(sequence, trinuc)
            percentage = calculate_percentage(sequence, trinuc)
            trinuc_results.append((trinuc, count, percentage))

        # Sort by percentage (descending)
        trinuc_results.sort(key=lambda x: x[2], reverse=True)

        output += f"{'Trinucleotide':<15}{'Count':<10}{'Percentage':<15}\n"
        output += "-" * 40 + "\n"

        for trinuc, count, percentage in trinuc_results:
            output += f"{trinuc:<15}{count:<10}{percentage:>6.2f}%\n"

        self.trinuc_text.insert(1.0, output)
        self.trinuc_results = trinuc_results

    def generate_summary(self, sequence):
        
        output = "=" * 70 + "\n"
        output += "SUMMARY STATISTICS\n"
        output += "=" * 70 + "\n\n"

        output += f"Sequence Analyzed: {sequence}\n"
        output += f"Sequence Length: {len(sequence)} nucleotides\n\n"

        
        output += "Nucleotide Composition:\n"
        output += "-" * 40 + "\n"
        for nuc in ['A', 'C', 'G', 'T']:
            count = sequence.count(nuc)
            percentage = (count / len(sequence)) * 100
            output += f"  {nuc}: {count} ({percentage:.2f}%)\n"

        output += "\n" + "=" * 70 + "\n"
        output += "TOP 5 DINUCLEOTIDES\n"
        output += "=" * 70 + "\n\n"

        for i, (dinuc, count, percentage) in enumerate(self.dinuc_results[:5], 1):
            output += f"  {i}. {dinuc}: {percentage:.2f}% ({count} occurrences)\n"

        output += "\n" + "=" * 70 + "\n"
        output += "TOP 5 TRINUCLEOTIDES\n"
        output += "=" * 70 + "\n\n"

        for i, (trinuc, count, percentage) in enumerate(self.trinuc_results[:5], 1):
            output += f"  {i}. {trinuc}: {percentage:.2f}% ({count} occurrences)\n"

        output += "\n" + "=" * 70 + "\n"

        self.summary_text.insert(1.0, output)

def main():
    root = tk.Tk()
    app = NucleotideAnalyzerGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
