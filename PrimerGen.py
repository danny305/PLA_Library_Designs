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


from PLA.RandSeqGen import *

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import MeltingTemp

import logging
import re
import math
from datetime import datetime
import json



def genPrimerPairs_5Ext(primer_length=20, anneal_length=10, GC_low=40, GC_high=60):
    """These primer pairs are designed so the template DNA will form a hair pin with 10 nucleotides
    as a 5' extension. This requires that the last 10 nucleotides in the reverse primer must be the
    first 10 nucleotides in the forward primer. """

    print('Primers for 5\' extension half-asstemers')

    forwTemplate5_3 = GenOligoGC(primer_length,GC_low, GC_high)
    """re.match checks if the first 2 Nuc are GC in the forward and backwards direction"""
    while not (re.match("[GC]{2}",str(forwTemplate5_3)) and
               re.match("[GC]{2}", str(forwTemplate5_3[::-1])) and
               re.match("[GC]{2}", str(forwTemplate5_3[10:12]))):

        forwTemplate5_3 = GenOligoGC(primer_length,GC_low, GC_high)

    forwTemp3_5 = forwTemplate5_3[::-1]
    forwPrimer5_3 = forwTemp3_5.complement()
    print(f"Template   Seq 3\' - > 5\': {forwTemp3_5}")
    print(f"ForwPrimer Seq 5\' - > 3\': {forwPrimer5_3}")

    forwPrimer_f10 = forwPrimer5_3[:10]
    print(f"First 10 Nucleotides of forward primer: {forwPrimer_f10}")

    revPrimer_f10 = GenOligoGC(10,GC_low, GC_high)
    while not re.match("[GC]{2}",str(revPrimer_f10)):
        revPrimer_f10 = GenOligoGC(10,GC_low, GC_high)

    revPrimer5_3 = revPrimer_f10 + forwPrimer_f10

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
               re.match("[GC]{2}", str(forwTemplate5_3[8:10]))):

        forwTemplate5_3 = GenOligoGC(primer_length,GC_low, GC_high)

    forwTemp3_5 = forwTemplate5_3[::-1]
    forwPrimer5_3 = forwTemp3_5.complement()
    print(f"Template   Seq 3\' - > 5\': {forwTemp3_5}")
    print(f"ForwPrimer Seq 5\' - > 3\': {forwPrimer5_3}")

    forwPrimer_L10 = forwPrimer5_3[10:]
    print(f"Last 10 Nucleotides of forward primer: {forwPrimer_L10}")

    revPrimer_L10 = GenOligoGC(10,GC_low, GC_high)
    while not re.match("[GC]{2}",str(revPrimer_L10[::-1])):
        revPrimer_L10 = GenOligoGC(10,GC_low, GC_high)

    """First 10 Nuc of rev primer must be identical to last 10 Nuc of forward Primer"""
    revPrimer5_3 = forwPrimer_L10 + revPrimer_L10

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

    fprimer_MT = MeltingTemp.Tm_GC(fprimer, Na=50, Mg=3, dNTPs=0.8)
    rprimer_MT = MeltingTemp.Tm_GC(rprimer, Na=50, Mg=3, dNTPs=0.8)

    fprimer_MT_NN = MeltingTemp.Tm_NN(fprimer, Na=50, Mg=3, dNTPs=0.8)
    rprimer_MT_NN = MeltingTemp.Tm_NN(fprimer, Na=50, Mg=3, dNTPs=0.8)


    print (f"forw primer: {fprimer}\nforw primer MT: {fprimer_MT} {fprimer_MT_NN} \n"
           f"rev  primer: {rprimer}\nrev primer MT : {rprimer_MT} {rprimer_MT_NN} \n")

    """Filters for primers that meet the MT standards"""
    if math.fabs(fprimer_MT - rprimer_MT) <= 3 and\
                max(fprimer_MT,rprimer_MT) <= 64 and\
                min(fprimer_MT, rprimer_MT) >= 60:

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






def main():
    pool = gen5_3PrimerPairPools()
    #convMastPool2File(pool, filename='testing', fext='txt')
    print(len(pool.keys()))
    for key, dict in pool.items():
        for key, value in dict.items():
            print(key, value)

main()


"""But in breaks and and check each function in between."""