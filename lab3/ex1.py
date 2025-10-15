#The melting temp is the temp at which one half of a particular DNA will disociate and become a single strand of DNA primerland and sequences are of critical importance in designing the parameters of a successful amplification. The tm of a nucleic acid duplex increases both with its length and with increasing GC content. tm=4(G+C)+2(A+T). Implement an app that computes the tm of a DNA sequence by using one of the formula or both of them. Input = a string of DNA, Output = temp in Celsius.
import math

def calc_temp(dna_seq):
    dna_seq = dna_seq.upper()
    a_count = dna_seq.count('A')
    t_count = dna_seq.count('T')
    c_count = dna_seq.count('C')
    g_count = dna_seq.count('G')
    Na = 0.001
    tm = 2 * (a_count + t_count) + 4 * (c_count + g_count)
    tm1 = 81.5 + 16.6 * (math.log(10)*Na) + 0.41 * ((c_count + g_count) / len(dna_seq) * 100) - (600 / len(dna_seq))
    return tm, tm1

dna_sequence = input("Enter a DNA sequence: ")
melting_temp = calc_temp(dna_sequence)
print(f"The melting temperature (Tm) of the DNA sequence is: {melting_temp} °C")

