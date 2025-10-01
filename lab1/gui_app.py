import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
from collections import Counter

class SequenceAnalyzerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Sequence Analyzer - DNA/RNA/Protein")
        self.root.geometry("800x600")

        
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill='both', expand=True, padx=10, pady=10)

        
        self.create_sequence_tab()

        
        self.create_fasta_tab()

    def create_sequence_tab(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Sequence Analysis")

        
        title = ttk.Label(frame, text="DNA/RNA/Protein Sequence Analyzer", font=('Arial', 14, 'bold'))
        title.pack(pady=10)

        
        input_frame = ttk.LabelFrame(frame, text="Enter Sequence", padding=10)
        input_frame.pack(fill='x', padx=20, pady=10)

        self.seq_input = tk.Text(input_frame, height=5, width=70)
        self.seq_input.pack(pady=5)
        self.seq_input.insert('1.0', 'ACGGGCATATGCGC')

       
        btn_frame = ttk.Frame(frame)
        btn_frame.pack(pady=10)

        analyze_btn = ttk.Button(btn_frame, text="Analyze Sequence", command=self.analyze_sequence)
        analyze_btn.pack(side='left', padx=5)

        clear_btn = ttk.Button(btn_frame, text="Clear", command=self.clear_sequence)
        clear_btn.pack(side='left', padx=5)

        
        results_frame = ttk.LabelFrame(frame, text="Results", padding=10)
        results_frame.pack(fill='both', expand=True, padx=20, pady=10)

        self.seq_results = scrolledtext.ScrolledText(results_frame, height=15, width=70)
        self.seq_results.pack(fill='both', expand=True)

    def create_fasta_tab(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="FASTA File Analysis")

        title = ttk.Label(frame, text="FASTA File Analyzer", font=('Arial', 14, 'bold'))
        title.pack(pady=10)

       
        file_frame = ttk.LabelFrame(frame, text="Select FASTA File", padding=10)
        file_frame.pack(fill='x', padx=20, pady=10)

        self.file_path_var = tk.StringVar()
        file_entry = ttk.Entry(file_frame, textvariable=self.file_path_var, width=60)
        file_entry.pack(side='left', padx=5)

        browse_btn = ttk.Button(file_frame, text="Browse", command=self.browse_file)
        browse_btn.pack(side='left', padx=5)

        analyze_file_btn = ttk.Button(file_frame, text="Analyze", command=self.analyze_fasta)
        analyze_file_btn.pack(side='left', padx=5)

        
        results_frame = ttk.LabelFrame(frame, text="Results", padding=10)
        results_frame.pack(fill='both', expand=True, padx=20, pady=10)

        self.fasta_results = scrolledtext.ScrolledText(results_frame, height=20, width=70)
        self.fasta_results.pack(fill='both', expand=True)

    def find_alphabet(self, sequence):
        return set(sequence)

    def analyze_composition(self, sequence):
        if not sequence:
            return {}
        alphabet = self.find_alphabet(sequence)
        percentages = {base: (sequence.count(base) / len(sequence)) * 100 for base in alphabet}
        return percentages

    def analyze_sequence(self):
        sequence = self.seq_input.get('1.0', 'end-1c').strip().upper()

        if not sequence:
            messagebox.showwarning("Warning", "Please enter a sequence!")
            return

        self.seq_results.delete('1.0', tk.END)

       
        self.seq_results.insert(tk.END, f"Sequence: {sequence}\n")
        self.seq_results.insert(tk.END, f"Length: {len(sequence)} characters\n\n")

        alphabet = self.find_alphabet(sequence)
        self.seq_results.insert(tk.END, f"Alphabet: {sorted(alphabet)}\n\n")

        composition = self.analyze_composition(sequence)
        self.seq_results.insert(tk.END, "Composition Analysis:\n")
        self.seq_results.insert(tk.END, "=" * 40 + "\n")

        for base in sorted(composition):
            self.seq_results.insert(tk.END, f"{base}: {composition[base]:.2f}%\n")

    def clear_sequence(self):
        self.seq_input.delete('1.0', tk.END)
        self.seq_results.delete('1.0', tk.END)

    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="Select FASTA File",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
        )
        if filename:
            self.file_path_var.set(filename)

    def read_fasta_file(self, file_path):
        sequences = []
        current_sequence = ""
        current_header = ""

        try:
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_sequence:
                            sequences.append((current_header, current_sequence))
                        current_header = line[1:]
                        current_sequence = ""
                    else:
                        current_sequence += line.upper()

                if current_sequence:
                    sequences.append((current_header, current_sequence))
        except FileNotFoundError:
            messagebox.showerror("Error", f"File '{file_path}' not found.")
            return None
        except Exception as e:
            messagebox.showerror("Error", f"Error reading file: {e}")
            return None

        return sequences

    def analyze_fasta(self):
        file_path = self.file_path_var.get()

        if not file_path:
            messagebox.showwarning("Warning", "Please select a FASTA file!")
            return

        self.fasta_results.delete('1.0', tk.END)

        sequences = self.read_fasta_file(file_path)

        if not sequences:
            self.fasta_results.insert(tk.END, "No sequences found in the file.\n")
            return

        self.fasta_results.insert(tk.END, f"Found {len(sequences)} sequence(s)\n\n")

        for i, (header, sequence) in enumerate(sequences, 1):
            self.fasta_results.insert(tk.END, "=" * 60 + "\n")
            self.fasta_results.insert(tk.END, f"Sequence {i}: {header}\n")
            self.fasta_results.insert(tk.END, "=" * 60 + "\n")
            self.fasta_results.insert(tk.END, f"Length: {len(sequence)} characters\n")

            alphabet = self.find_alphabet(sequence)
            self.fasta_results.insert(tk.END, f"Alphabet: {sorted(alphabet)}\n\n")

            composition = self.analyze_composition(sequence)
            self.fasta_results.insert(tk.END, "Composition Analysis:\n")
            self.fasta_results.insert(tk.END, "-" * 40 + "\n")

            for base in sorted(composition):
                self.fasta_results.insert(tk.END, f"{base}: {composition[base]:.2f}%\n")

            preview = sequence[:50] + ('...' if len(sequence) > 50 else '')
            self.fasta_results.insert(tk.END, f"\nSequence preview: {preview}\n\n")

if __name__ == "__main__":
    root = tk.Tk()
    app = SequenceAnalyzerGUI(root)
    root.mainloop()
