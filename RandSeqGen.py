from random import choice, seed
from math import pow, log
from pprint import pprint
from time import time

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import MeltingTemp as MT

import re



# Calculate the number of nucleotides needed for a diverse library
def oligoLenCalc(magnitude=14):
	oligo_len = log(pow(10, magnitude), 4)
	return oligo_len

# Generate random sequence of DNA
def randSeqGen(length=25):

    seq = [choice(['A', 'T', 'C', 'G']) for _ in range(length)]
    seq = Seq(''.join(seq), generic_dna)
    if re.search("[G]{4,100}|[C]{4,100}|[A]{4,100}|[T]{4,100}]", str(seq)):
        #print(seq)
        #print('This sequence has at least 4 consecutive same nucleotides.')
        seq = randSeqGen(length)

    return seq


def chkGC_content(sequence=Seq('GAGACCTT',generic_dna)):
    sequence = sequence.upper()
    seq_len = len(sequence)
    num_G, num_C = sequence.count('G'), sequence.count('C')
    GC_content = round((num_C + num_G) / seq_len * 100, 1)
    return GC_content

def GenOligoGC(length=25, GC_low=40, GC_high=60):
    generate = True
    seed(time())
    while generate == True:
        seq = randSeqGen(length)
        GC_content = chkGC_content(seq)
        if GC_content >= GC_low and GC_content <= GC_high:
            return seq

def printNSeq(num_seq=20,seq_len=20,GC_low_cutoff=40,GC_high_cutoff=60):
    uniq_seq = list()
    while len(uniq_seq) < num_seq:
        seq = randSeqGen(seq_len)
        GCcontent = chkGC_content(seq)
        if GCcontent >= GC_low_cutoff and GCcontent <= GC_high_cutoff:
            uniq_seq.append((str(seq), GCcontent, round(MT.Tm_GC(seq, Na=50, Mg=3, dNTPs=0.8),3)))
        else:
            continue
    pprint(uniq_seq)

def main():
    num_uniq_seq = int(input('How many unique sequences do you want to generate:'))
    seq_len = int(input('What is the desired sequence length:'))
    GC_low_cutoff = round(float(input('GC content lower cut-off(up to 2 decimal):')), 2)
    GC_high_cutoff = round(float(input('GC content upper cut-off(up to 2 decimal):')), 2)

    uniq_seq = list()
    while len(uniq_seq) < num_uniq_seq:
        seq = randSeqGen(seq_len)
        GCcontent = chkGC_content(seq)
        if GCcontent >= GC_low_cutoff and GCcontent <= GC_high_cutoff:
            uniq_seq.append((str(seq),GCcontent))
        else:
            continue
    pprint(uniq_seq)

if __name__=="__main__":
    #main()
    printNSeq()
