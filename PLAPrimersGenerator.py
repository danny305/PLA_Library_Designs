"""PLA Primer Generator 5' and 3' extension Sequences"""
"""
November 2nd, 2018

Notes for Good Primer Design..
 - Length of 18-24 bases                                        check
 - 40-60% G/C content                                           check                                 
 - Start and end with 1-   2 G/C pairs                          check
 - Melting temperature (Tm) of 50-60°C                          check
 - Primer pairs should have a Tm within 5°C of each other
 - Primer pairs should not have complementary regions           This partially applies to us. 
"""

"""
Folding primer complementation
"""


"""IDT Primer Recommendations

IDT recommends that you aim for PCR primers between 18 and 30 bases; however, the most 
important considerations for primer design should be their Tm value and specificity. 
Primers should also be free of strong secondary structures and self-complementarity. 
Design your PCR primers to conform to the following guidelines:

Melting temperature (Tm): The optimal melting temperature of the primers is 60–64°C, with 
an ideal temperature of 62°C, which is based on typical cycling and reaction conditions 
and the optimum temperature for PCR enzyme function. Ideally, the melting temperatures of 
the 2 primers should not differ by more than 2°C in order for both primers to bind 
simultaneously and efficiently amplify the product.

Annealing temperature (Ta): The annealing temperature chosen for PCR relies directly on 
length and composition of the primers. This temperature should be no more than 5°C below 
the Tm of your primers. One consequence of having Ta too low is that one or both primers 
will anneal to sequences other than the intended target because internal single-base 
mismatches or partial annealing may be tolerated. This can lead to nonspecific PCR 
amplification and will consequently reduce the yield of the desired product. Conversely, 
if Ta is too high, reaction efficiency may be reduced because the likelihood of primer 
annealing is reduced significantly. Optimal annealing temperatures will result in the 
highest product yield with the correct amplicon.

GC content: Design your assay so that the GC content is 35–65%, with an ideal content of 
50%, which allows complexity while still maintaining a unique sequence. Primer sequences 
should not contain regions of 4 or more consecutive G residues.

"""

from random import choice, seed
from math import pow, log
from pprint import pprint
from time import time

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import MeltingTemp as MT
from Bio.SeqUtils import GC as GCpercent

import logging
import re
import math
from datetime import datetime
import json



# Calculate the number of nucleotides needed for a diverse library
def oligoLenCalc(magnitude=14):
	oligo_len = log(pow(10, magnitude), 4)
	return oligo_len

# Generate random sequence of DNA
def randSeqGen(length=25):
    """Generates and returns a DNA Sequence object of the desired nucleotide length.
    Where the maximum number of adjacent same nucleotides is 3"""
    seq_approved = False

    while not seq_approved:
        seq = [choice(['A', 'T', 'C', 'G']) for _ in range(length)]
        seq = Seq(''.join(seq), generic_dna)
        if re.search("[G]{4,100}|[C]{4,100}|[A]{4,100}|[T]{4,100}]", str(seq)):
            #print(f"The sequence: {seq} is invalid")
            continue

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
        GC_content = GCpercent(seq)
        if GC_content >= GC_low and GC_content <= GC_high:
            return seq



def printNSeq(num_seq=20,seq_len=20,GC_low_cutoff=40,GC_high_cutoff=60):
    uniq_seq = list()
    while len(uniq_seq) < num_seq:
        seq = randSeqGen(seq_len)
        GCcontent = GCpercent(seq)
        if GCcontent >= GC_low_cutoff and GCcontent <= GC_high_cutoff:
            uniq_seq.append((str(seq), GCcontent, round(MT.Tm_GC(seq, Na=50, Mg=3, dNTPs=0.8),3)))
        else:
            continue
    pprint(uniq_seq)



def genPrimerPairs_5Ext(primer_length=20, anneal_length=10, GC_low=40, GC_high=60):
    """These primer pairs are designed so the template DNA will form a hair pin with 10 nucleotides
    as a 5' extension. This requires that the last 10 nucleotides in the reverse primer must be the
    first 10 nucleotides in the forward primer. """

    print('Primers for 5\' extension half-asstemers')

    forwTemplate5_3 = GenOligoGC(primer_length,GC_low, GC_high)
    """re.match checks if the first 2 Nuc are GC in the forward and backwards direction"""
    while not (re.match("[GC]{2}",str(forwTemplate5_3)) and
               re.match("[GC]{2}", str(forwTemplate5_3[::-1])) and
               re.match("[GC]{2}", str(forwTemplate5_3[anneal_length:anneal_length+2]))):

        forwTemplate5_3 = GenOligoGC(primer_length,GC_low, GC_high)

    forwTemp3_5 = forwTemplate5_3[::-1]
    forwPrimer5_3 = forwTemp3_5.complement()
    print(f"Template   Seq 3\' - > 5\': {forwTemp3_5}")
    print(f"ForwPrimer Seq 5\' - > 3\': {forwPrimer5_3}")

    forwPrimer_fHalf = forwPrimer5_3[:anneal_length]
    print(f"First 10 Nucleotides of forward primer: {forwPrimer_fHalf}")

    revPrimer_fHalf = GenOligoGC(anneal_length,GC_low, GC_high)
    while not re.match("[GC]{2}",str(revPrimer_fHalf)):
        revPrimer_fHalf = GenOligoGC(anneal_length,GC_low, GC_high)

    revPrimer5_3 = revPrimer_fHalf + forwPrimer_fHalf

    print(f"RevPrimer  Seq 5\' - > 3\': {revPrimer5_3}")

    return forwPrimer5_3, revPrimer5_3



def genPrimerPairs_3Ext(primer_length=20, anneal_length=10, GC_low=40, GC_high=60):
    """These primer pairs are designed so the template DNA will form a hair pin with 10 nucleotides
    as a 3' extension. This requires that the first 10 nucleotides in the reverse primer must be
    identical to the last 10 nucleotides in the forward primer. """

    print('Primers for 3\' extension half-asstemers')


    forwTemplate5_3 = GenOligoGC(primer_length,GC_low, GC_high)
    """re.match checks if the first 2 Nuc are GC in the forward and backwards direction"""
    while not (re.match("[GC]{2}",str(forwTemplate5_3)) and
               re.match("[GC]{2}", str(forwTemplate5_3[::-1])) and
               re.match("[GC]{2}", str(forwTemplate5_3[anneal_length-2:anneal_length]))):

        forwTemplate5_3 = GenOligoGC(primer_length,GC_low, GC_high)

    forwTemp3_5 = forwTemplate5_3[::-1]
    forwPrimer5_3 = forwTemp3_5.complement()
    print(f"Template   Seq 3\' - > 5\': {forwTemp3_5}")
    print(f"ForwPrimer Seq 5\' - > 3\': {forwPrimer5_3}")

    forwPrimer_LHalf = forwPrimer5_3[anneal_length:]
    print(f"Last 10 Nucleotides of forward primer: {forwPrimer_LHalf}")

    revPrimer_LHalf = GenOligoGC(anneal_length,GC_low, GC_high)
    while not re.match("[GC]{2}",str(revPrimer_LHalf[::-1])):
        revPrimer_LHalf = GenOligoGC(anneal_length,GC_low, GC_high)

    """First half Nucleotides of rev primer must be identical to last half Nucleotides of forward Primer"""
    revPrimer5_3 = forwPrimer_LHalf + revPrimer_LHalf

    print(f"RevPrimer  Seq 5\' - > 3\': {revPrimer5_3}")

    return forwPrimer5_3, revPrimer5_3



def evalPrimerPairMT(fprimer, rprimer, ret_mt=False):
    """This will check the melting temperature

    The optimal melting temperature of the primers is 60–64°C, with
    an ideal temperature of 62°C, which is based on typical cycling and reaction conditions
    and the optimum temperature for PCR enzyme function. Ideally, the melting temperatures of
    the 2 primers should not differ by more than 2°C in order for both primers to bind
    simultaneously and efficiently amplify the product.
    PCR parameters used are from IDT: Oligo 0.2 uM Na 50 mM, Mg 3 mM, dNTPs 0.8 mM
    :param ret_mt: """

    fprimer_MT = MT.Tm_GC(fprimer, Na=50, Mg=3, dNTPs=0.8)
    rprimer_MT = MT.Tm_GC(rprimer, Na=50, Mg=3, dNTPs=0.8)

    fprimer_MT_NN = MT.Tm_NN(fprimer, Na=50, Mg=3, dNTPs=0.8)
    rprimer_MT_NN = MT.Tm_NN(fprimer, Na=50, Mg=3, dNTPs=0.8)


    print (f"forw primer: {fprimer}\nforw primer MT: {fprimer_MT} {fprimer_MT_NN} \n"
           f"rev  primer: {rprimer}\nrev primer MT : {rprimer_MT} {rprimer_MT_NN} \n")

    """Filters for primers that meet the MT standards"""
    if math.fabs(fprimer_MT - rprimer_MT) <= 3.0 and\
                max(fprimer_MT,rprimer_MT) <= 64.0 and\
                min(fprimer_MT, rprimer_MT) >= 59.0:

        print("MT of primer pair passed.\n")

        if ret_mt == False:
            return True
        else:
            return fprimer_MT, rprimer_MT
    else:
        print("MT for the primer pairs did not meet standards\n"); return False



def retPrimerPairMT(fprimer,rprimer):
    forwMT, revMT = [round(temp,2) for temp in evalPrimerPairMT(
                                                fprimer, rprimer, ret_mt=True)]

    return forwMT, revMT



def genCertPrimerPairs(length=20,*,GC_low=40, GC_high=60,ext=5, ret_str=True):

    chk_primer_mt = False
    while not chk_primer_mt:

        if ext == 5:
            fprimer, rprimer = genPrimerPairs_5Ext(primer_length=length, GC_low=GC_low, GC_high=GC_high)
        elif  ext == 3:
            fprimer, rprimer = genPrimerPairs_3Ext(primer_length=length, GC_low=GC_low, GC_high=GC_high)

        chk_primer_mt = evalPrimerPairMT(fprimer, rprimer)

    if ret_str == False:
        return fprimer,rprimer
    else:
        return str(fprimer), str(rprimer)



def genPrimerPairPool(pool_size=8,length=20,*,GC_low=40, GC_high=60,ext=5, ret_str=True):
    """Generate a pool of primer pairs that pass the melt temp criteria defined in
    genCertPrimerPairs(). Default is the 5'extension half-asstemers primers."""

    pool = dict()

    while len(pool.keys()) <= pool_size:
        pool['meta_data'] = {'pool_size': pool_size, 'primer_length': length,
                             "GC_low": GC_low, "GC_high": GC_high, "extension_end": ext
                             }
        primerKeyName = "primers " + str(len(pool.keys()))
        primers = dict()
        primers['forw_primer'], primers['rev_primer'] = \
                                    genCertPrimerPairs(length=length,GC_low=GC_low,
                                                       GC_high=GC_high, ext=ext,
                                                       ret_str=ret_str)

        primers['forw_MT'], primers['rev_MT'] = retPrimerPairMT(primers['forw_primer'],
                                                                primers['rev_primer'])
        pool[primerKeyName] = primers


    return pool



def gen5_3PrimerPairPools(pool_size=8,length=20,*,GC_low=40, GC_high=60,ret_str=True):
    """Return primer pair pools for the 5' extension and the 3' extension."""

    masterPool = dict()
    pool_5Ext = genPrimerPairPool(pool_size=pool_size,length=length,GC_low=GC_low,
                                  GC_high=GC_high,ext=5,ret_str=ret_str)

    pool_3Ext = genPrimerPairPool(pool_size=pool_size, length=length, GC_low=GC_low,
                                  GC_high=GC_high, ext=3, ret_str=ret_str)

    masterPool["5\'_ext_pool"] = pool_5Ext; masterPool["3\'_ext_pool"] = pool_3Ext

    return masterPool



def convMastPool2File(mast_pool_dict,*,filename,fext='json'):
    now = datetime.now().strftime("%m-%d-%y_%H:%M")
    with open(f"{filename}-{now}.{fext}",'w+') as f:
        json.dump(mast_pool_dict, f, indent=4,sort_keys=True)
    print('created file')



"""Primers that can self-hybridize will be unavailable for hybridization to the template. 
Generally avoid primers that can form 4 or more consecutive bonds with itself, or 8 or more 
bonds total. Example of a marginally problematic primer:

This oligo forms a substantially stable dimer with itself, with four consecutive bonds at two places 
and a total of eight inter-strand bonds. Primers with 3' ends hybridizing even transiently will become 
extended due to polymerase action, thus ruining the primer and generating false bands. Be somewhat more 
stringent in avoiding 3' dimers.
"""


def pullAllPrimers(pool_dict=None):
    pass





def main():
    #gen5_3PrimerPairPools()
    pullAllPrimers()


if __name__=="__main__":
    main()