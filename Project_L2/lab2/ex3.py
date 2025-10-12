#Design an application using the AI which uses a GUI and allows the user to selecte a fasta file. The content of the file should be analyzed by using a sliding windoow of 30 positions. The content for each sliding window should be used in order to extract the relative freqeuncies of the symbols found in the alphabet of the sequence which is the content of the fasta file. The input will be the dna sequence of the fasta file and the output should be the values of the relative frequencies of each symbol in the alphabet of the sequence. Translate it as lines on a chart. Thus your chart should have 4 lines which reflect the values found over the sequence
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

VALID = ("A", "C", "G", "T")

def read_fasta_first(path: str) -> str:
    seq = []
    with open(path, "r", encoding="utf-8") as f:
        started = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if started and seq:
                    break
                started = True
                continue
            if started:
                seq.append(line)
    return "".join(seq).upper()

def rolling_relative_freqs(seq: str, k: int):
    n = len(seq)
    if k <= 0 or n < k:
        return [], {b: [] for b in VALID}

    def contrib(ch):
        return (
            1 if ch == "A" else 0,
            1 if ch == "C" else 0,
            1 if ch == "G" else 0,
            1 if ch == "T" else 0,
            1 if ch in VALID else 0,  
        )

    a_sum = c_sum = g_sum = t_sum = denom = 0
    for i in range(k):
        a, c, g, t, v = contrib(seq[i])
        a_sum += a; c_sum += c; g_sum += g; t_sum += t; denom += v

    xs = []
    freqs = {b: [] for b in VALID}

    def record(center_idx, a_sum, c_sum, g_sum, t_sum, denom):
        xs.append(center_idx)
        if denom == 0:
            for b in VALID: freqs[b].append(0.0)
        else:
            freqs["A"].append(a_sum / denom)
            freqs["C"].append(c_sum / denom)
            freqs["G"].append(g_sum / denom)
            freqs["T"].append(t_sum / denom)

    half = (k - 1) / 2.0
    record(1 + half, a_sum, c_sum, g_sum, t_sum, denom)

    for start in range(1, n - k + 1):
        out_ch = seq[start - 1]
        in_ch  = seq[start + k - 1]

        ao, co, go, to, vo = contrib(out_ch)
        a_sum -= ao; c_sum -= co; g_sum -= go; t_sum -= to; denom -= vo

        ai, ci, gi, ti, vi = contrib(in_ch)
        a_sum += ai; c_sum += ci; g_sum += gi; t_sum += ti; denom += vi

        record(start + 1 + half, a_sum, c_sum, g_sum, t_sum, denom)

    return xs, freqs

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Sliding-Window Relative Frequencies (DNA)")
        self.geometry("1000x680")
        self.seq = ""
        self.filepath = ""

        top = ttk.Frame(self, padding=10)
        top.pack(side=tk.TOP, fill=tk.X)

        ttk.Button(top, text="Open FASTA…", command=self.open_file).pack(side=tk.LEFT)
        self.file_var = tk.StringVar(value="No file loaded")
        ttk.Label(top, textvariable=self.file_var, width=50, anchor="w").pack(side=tk.LEFT, padx=10)

        ttk.Label(top, text="Window size:").pack(side=tk.LEFT, padx=(20,5))
        self.k_var = tk.StringVar(value="30")
        ttk.Entry(top, textvariable=self.k_var, width=6).pack(side=tk.LEFT)

        ttk.Button(top, text="Analyze & Plot", command=self.run).pack(side=tk.RIGHT)

        info = ttk.Frame(self, padding=(10,0))
        info.pack(side=tk.TOP, fill=tk.X)
        self.info_var = tk.StringVar(value="Load a FASTA file to begin.")
        ttk.Label(info, textvariable=self.info_var, anchor="w").pack(side=tk.LEFT)

        self.fig = Figure(figsize=(10,4.8), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel("Window center (1-based position)")
        self.ax.set_ylabel("Relative frequency")
        self.ax.set_ylim(0, 1)
        self.ax.grid(True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def open_file(self):
        path = filedialog.askopenfilename(
            title="Choose FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffa *.txt"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            self.seq = read_fasta_first(path)
        except Exception as e:
            messagebox.showerror("Error", f"Could not read FASTA: {e}")
            return
        if not self.seq:
            messagebox.showwarning("Empty", "No sequence found in the FASTA file.")
            return
        self.filepath = path
        self.file_var.set(os.path.basename(path))
        self.info_var.set(f"Sequence length: {len(self.seq):,} bases. Non-ACGT ignored in denominator.")

    def run(self):
        if not self.seq:
            messagebox.showinfo("No sequence", "Please open a FASTA file first.")
            return
        try:
            k = int(self.k_var.get())
            if k <= 0:
                raise ValueError
        except Exception:
            messagebox.showerror("Invalid window", "Window size must be a positive integer.")
            return

        xs, freqs = rolling_relative_freqs(self.seq, k)
        if not xs:
            messagebox.showinfo("Too short", f"Sequence is shorter than window size ({k}).")
            return

        self.ax.clear()
        self.ax.set_title(f"Sliding-window (k={k}) relative frequencies")
        self.ax.set_xlabel("Window center (1-based position)")
        self.ax.set_ylabel("Relative frequency")
        self.ax.set_ylim(0, 1)
        self.ax.grid(True)

        self.ax.plot(xs, freqs["A"], label="A")
        self.ax.plot(xs, freqs["C"], label="C")
        self.ax.plot(xs, freqs["G"], label="G")
        self.ax.plot(xs, freqs["T"], label="T")
        self.ax.legend(loc="upper right")

        self.canvas.draw()

if __name__ == "__main__":
    App().mainloop()

