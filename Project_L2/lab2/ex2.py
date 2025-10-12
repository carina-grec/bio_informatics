#Find in sequence S only the dinucleo and trinucleo that exist, without the use of the bruteforce engine. In order to achive the results, one must check these combinations starting from the beginning of the sequence until the end of the sequence.
S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

S = S.upper()

din = set()
trin = set()

for i in range(len(S) - 1):
    di = S[i:i+2]
    if len(di) == 2 and all(base in "ACGT" for base in di):
        din.add(di)
for i in range(len(S) - 2):
    tri = S[i:i+3]
    if len(tri) == 3 and all(base in "ACGT" for base in tri):
        trin.add(tri)

print("Dinucleotides found in sequence:")
print(sorted(din))

print("\nTrinucleotides found in sequence:")
print(sorted(trin))

