import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt

def parse_fasta(file_path):
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    sequence = ''
    for line in lines:
        if not line.startswith('>'):
            sequence += line.strip()
    return sequence.upper()

def compute_frequencies(sequence, window_size=30):
    
    freq_A = []
    freq_C = []
    freq_G = []
    freq_T = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        a = window.count('A') / window_size
        c = window.count('C') / window_size
        g = window.count('G') / window_size
        t = window.count('T') / window_size
        freq_A.append(a)
        freq_C.append(c)
        freq_G.append(g)
        freq_T.append(t)
    return freq_A, freq_C, freq_G, freq_T

class FrequencyGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Nucleotide Frequency Analysis")
        self.file_path = None

      
        self.select_btn = tk.Button(root, text="Select FASTA File", command=self.select_file)
        self.select_btn.pack(pady=10)

       
        self.analyze_btn = tk.Button(root, text="Analyze and Plot", command=self.analyze)
        self.analyze_btn.pack(pady=10)

    def select_file(self):
        self.file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])

    def analyze(self):
        if not self.file_path:
            messagebox.showerror("Error", "Please select a FASTA file first.")
            return

        try:
            sequence = parse_fasta(self.file_path)
            if not sequence:
                messagebox.showerror("Error", "No sequence found in the file.")
                return

            freq_A, freq_C, freq_G, freq_T = compute_frequencies(sequence)

            
            plt.figure(figsize=(10, 6))
            plt.plot(freq_A, label='A', color='red')
            plt.plot(freq_C, label='C', color='blue')
            plt.plot(freq_G, label='G', color='green')
            plt.plot(freq_T, label='T', color='orange')
            plt.xlabel('Window Position')
            plt.ylabel('Relative Frequency')
            plt.title('Nucleotide Frequencies in Sliding Windows (Size 30)')
            plt.legend()
            plt.grid(True)
            plt.show()

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = FrequencyGUI(root)
    root.mainloop()
