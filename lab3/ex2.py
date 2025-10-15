#Create an app that uses the sliding window method in order to read the melting temp (tm) over the sequence s. Use a sliding window of 8 positions and choose the fasta file as input.  
#Do the plotting graph for the sliding window.

import math
import matplotlib.pyplot as plt

def sw_tm(dna_seq, window_size=8):
    tm_val = []
    for i in range(len(dna_seq) - window_size + 1):
        win = dna_seq[i:i + window_size]
        a_count = win.count('A')
        t_count = win.count('T')
        c_count = win.count('C')
        g_count = win.count('G')
        tm = 2 * (a_count + t_count) + 4 * (c_count + g_count)
        tm_val.append((i, tm))
    return tm_val

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
        return seq
    
def plot_tm(tm_values):
    positions, tm = zip(*tm_values)
    plt.plot(positions, tm)
    plt.xlabel('Position')
    plt.ylabel('Melting Temperature (Tm)')
    plt.title('Sliding Window Melting Temperature')
    plt.show()
fasta_file = input("Enter the path to the FASTA file: ")
sequence = read_fasta(fasta_file)
tm_values = sw_tm(sequence)
plot_tm(tm_values)




