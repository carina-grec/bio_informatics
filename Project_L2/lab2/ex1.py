#Find the percentage for all the dinucleotide and trinucleotide combinations for the sequence:
#1. Build a brute force engine to generate all dinucleotide and trinucleotide combinations.
#2. For each combination, find out the percentage inside the S sequence.
#3. Show the percentage for each combination in the output of your implementation.
from itertools import product

# your sequence
S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
bases = ['A', 'C', 'G', 'T']

dinucs  = [''.join(p) for p in product(bases, repeat=2)]
trinucs = [''.join(p) for p in product(bases, repeat=3)]

def percentage(seq, combo):
    seq = seq.upper()
    count = seq.count(combo)
    possible = len(seq) - len(combo) + 1
    return (count / possible * 100) if possible > 0 else 0

print("Dinucleotide percentages:")
for d in dinucs:
    print(f"{d}: {percentage(S, d):.2f}%")

print("\nTrinucleotide percentages:")
for t in trinucs:
    print(f"{t}: {percentage(S, t):.2f}%")


