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
from collections import namedtuple
from difflib import SequenceMatcher
from itertools import combinations

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
        pool[f'meta_data_{ext}'] = {'pool_size': pool_size, 'primer_length': length,
                             "GC_low": GC_low, "GC_high": GC_high, "extension_end": ext
                             }
        primerNum = str(len(pool.keys()))
        primerKeyName = f"primers {ext}.{primerNum}"
        primers = dict()
        primers[f'forw_primer_{ext}.{primerNum}'], primers[f'rev_primer_{ext}.{primerNum}'] = \
            genCertPrimerPairs(length=length,GC_low=GC_low, GC_high=GC_high, ext=ext, ret_str=ret_str)

        primers[f'forw_MT_{ext}.{primerNum}'], primers[f'rev_MT_{ext}.{primerNum}'] = \
            retPrimerPairMT(primers[f'forw_primer_{ext}.{primerNum}'], primers[f'rev_primer_{ext}.{primerNum}'])

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
    if pool_dict == None:
        pool_dict = gen5_3PrimerPairPools(50)

    dict1, dict2 = [ dict for dict in pool_dict.values()]
    all_primers = {**dict1,**dict2}
    keys_2_pop = [key for key in all_primers.keys() if re.match('meta.*',key)]
    for meta_key in keys_2_pop:
        all_primers.pop(meta_key)

    all_Seq = list()
    PrimerPair = namedtuple('PrimerPair', ['extension','forw_primer','rev_primer', 'forw_MT','rev_MT'])

    for primer_pair in all_primers.keys():
        primer = all_primers[primer_pair]
        if re.search('5\.',list(primer.keys())[0]):
            ext='5_prime'
        elif re.search('3\.',list(primer.keys())[0]):
            ext='3_prime'

        all_Seq.append(PrimerPair(ext,*primer.values()))

    return all_Seq


def chkSelfDimerization(all_seq):

    filtered_seq = list()
    for (index, (_, forwP, revP, *MTs)) in enumerate(all_seq):
        print(index)
        forwP,revP = Seq(forwP,generic_dna), Seq(revP,generic_dna)
        forwP_Rev = forwP[::-1]
        forwP_Com = forwP.complement()
        match_ForwP = SequenceMatcher(a=forwP_Rev,b=forwP_Com).find_longest_match(0,len(forwP_Rev),0,len(forwP_Com))
        match_forwP_block = SequenceMatcher(a=forwP_Rev, b=forwP_Com).get_matching_blocks()
        print(match_ForwP)
        print(match_forwP_block)
        print(forwP_Rev[match_ForwP.a:match_ForwP.a+match_ForwP.size], forwP_Com[match_ForwP.b:match_ForwP.b+match_ForwP.size],sep='\n')
        print("Forw_Rev: ",forwP_Rev, "Complement: ",forwP.complement(),sep='\n',end='\n\n')

        revP_Rev = revP[::-1]
        revP_Com = revP.complement()
        match_RevP = SequenceMatcher(a=revP_Rev, b=revP_Com).find_longest_match(0, len(revP_Rev), 0, len(revP_Com))
        match_RevP_block = SequenceMatcher(a=revP_Rev, b=revP_Com).get_matching_blocks()
        print("Reverse Primer")
        print(match_RevP)
        print(match_RevP_block)
        print(revP_Rev[match_RevP.a:match_RevP.a + match_RevP.size], revP_Com[match_RevP.b:match_RevP.b + match_RevP.size], sep='\n')
        print("RevP_Rev: ", revP_Rev, "Complement: ", revP.complement(), sep='\n', end='\n\n')

        if match_ForwP.size > 3 or match_RevP.size > 3:
            continue
        else:
            print(f'    Adding index: {index}',end='\n\n')
            filtered_seq.append(all_seq[index])

    print(filtered_seq)
    print(len(filtered_seq),end='\n\n')
    return filtered_seq


def seqCombinations(all_seq):
    print(250*"_")
    all_combs = combinations(all_seq,2)
    for (index,(pair1, pair2)) in enumerate(all_combs):
        _, pair1_forwP, pair1_revP, *MTs = pair1
        _, pair2_forwP, pair2_revP, *MTs = pair2

        pair1_forwP, pair1_revP = Seq(pair1_forwP, generic_dna), Seq(pair1_revP, generic_dna)
        pair2_forwP, pair2_revP = Seq(pair2_forwP, generic_dna), Seq(pair2_revP, generic_dna)

        pair1_forwP_Com, pair1_revP_Com = pair1_forwP.complement(), pair1_revP.complement()
        pair2_forwP_Rev, pair2_revP_Rev = pair2_forwP[::-1], pair2_revP[::-1]


        """ Pair1_Forward Primer Being compared to Pair2_Forward Primer"""

        match_both_ForwP = SequenceMatcher(a=pair1_forwP_Com, b=pair2_forwP_Rev).find_longest_match(0, len(pair2_forwP_Rev), 0, len(pair1_forwP_Com))
        match_both_forwP_block = SequenceMatcher(a=pair1_forwP_Com, b=pair2_forwP_Rev).get_matching_blocks()

        print(index)
        print(match_both_ForwP)
        print(match_both_forwP_block)
        print(pair2_forwP_Rev[match_both_ForwP.b:match_both_ForwP.b + match_both_ForwP.size],
              pair1_forwP_Com[match_both_ForwP.a:match_both_ForwP.a + match_both_ForwP.size], sep='\n')
        print("Pair2_ForwP_Rev: ", pair2_forwP_Rev, "Pair1_ForwP_Complement: ", pair1_forwP_Com, sep='\n', end='\n\n')




        """ Pair1_Forward Primer Being compared to Pair2_Reverse Primer"""

        match_P1ForwP_P2Rev = SequenceMatcher(a=pair1_forwP_Com, b=pair2_revP_Rev).find_longest_match(0, len(
            pair2_revP_Rev), 0, len(pair1_forwP_Com))
        match_P1ForwP_P2Rev_block = SequenceMatcher(a=pair1_forwP_Com, b=pair2_revP_Rev).get_matching_blocks()

        print(match_P1ForwP_P2Rev)
        print(match_P1ForwP_P2Rev_block)
        print(pair2_revP_Rev[match_P1ForwP_P2Rev.b:match_P1ForwP_P2Rev.b + match_P1ForwP_P2Rev.size],
              pair1_forwP_Com[match_P1ForwP_P2Rev.a:match_P1ForwP_P2Rev.a + match_P1ForwP_P2Rev.size], sep='\n')
        print("Pair2_RevP_Rev: ", pair2_revP_Rev, "Pair1_ForwP_Complement: ", pair1_forwP_Com, sep='\n', end='\n\n')





        """ Pair1_Reverse Primer Being compared to Pair2_Forward Primer"""

        match_P1RevP_P2ForwP = SequenceMatcher(a=pair1_revP_Com, b=pair2_forwP_Rev).find_longest_match(0, len(
            pair2_forwP_Rev), 0, len(pair1_revP_Com))
        match_P1RevP_P2Forw_block = SequenceMatcher(a=pair1_revP_Com, b=pair2_forwP_Rev).get_matching_blocks()

        print(match_P1RevP_P2ForwP)
        print(match_P1RevP_P2Forw_block)
        print(pair2_forwP_Rev[match_P1RevP_P2ForwP.b:match_P1RevP_P2ForwP.b + match_P1RevP_P2ForwP.size],
              pair1_revP_Com[match_P1RevP_P2ForwP.a:match_P1RevP_P2ForwP.a + match_P1RevP_P2ForwP.size], sep='\n')
        print("Pair2_forwP_Rev: ", pair2_forwP_Rev, "Pair1_RevP_Complement: ", pair1_revP_Com, sep='\n', end='\n\n')






        """ Pair1_Reverse Primer Being compared to Pair2_Reverse Primer"""

        match_P1RevP_P2RevP = SequenceMatcher(a=pair1_revP_Com, b=pair2_revP_Rev).find_longest_match(0, len(
            pair2_revP_Rev), 0, len(pair1_revP_Com))
        match_P1RevP_P2Rev_block = SequenceMatcher(a=pair1_revP_Com, b=pair2_revP_Rev).get_matching_blocks()

        print(match_P1RevP_P2RevP)
        print(match_P1RevP_P2Rev_block)
        print(pair2_revP_Rev[match_P1RevP_P2RevP.b:match_P1RevP_P2RevP.b + match_P1RevP_P2RevP.size],
              pair1_revP_Com[match_P1RevP_P2RevP.a:match_P1RevP_P2RevP.a + match_P1RevP_P2RevP.size], sep='\n')
        print("Pair2_RevP_Rev: ", pair2_revP_Rev, "Pair1_RevP_Complement: ", pair1_revP_Com, sep='\n', end='\n\n')

        """
        match_both_ForwP
        match_P1ForwP_P2Rev
        match_P1RevP_P2ForwP
        match_P1RevP_P2RevP
        """
        #ToDo These need to be compared and need to decide how to discard and keep the correct primers


        # if match_ForwP.size > 3 or match_RevP.size > 3:
        #     continue
        # else:
        #     print(f'    Adding index: {index}',end='\n\n')
        #     filtered_seq.append(all_seq[index])






        #print(pair1, pair2, sep="-----")









def main():
    # mast_dict = gen5_3PrimerPairPools()
    # convMastPool2File(mast_dict,filename="test", fext='txt')
    all_seq =pullAllPrimers()
    filt_seq = chkSelfDimerization(all_seq)
    seqCombinations(filt_seq)



if __name__=="__main__":
    main()