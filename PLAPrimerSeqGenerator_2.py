


from numpy.random import choice, seed
from math import pow, log
from pprint import pprint
from time import time
from collections import namedtuple
from difflib import SequenceMatcher
from itertools import combinations

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import MeltingTemp as MT
from Bio.SeqUtils import GC

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

        if re.search("[G]{4,100}|[C]{4,100}|[A]{4,100}|[T]{4,100}", str(seq)):
            #print(f"The sequence: {seq} is invalid")
            continue

        if len(re.findall('C{2,3}|G{2,3}|A{2,3}|T{2,3}',str(seq))) > 4 or \
            len(re.findall('C{3}|G{3}|A{3}|T{3}', str(seq))) > 2:
            #print(f'sequence not random enough: {seq}')
            continue

        print(re.findall('[C]{2,3}|[G]{2,3}|[A]{2,3}|[T]{2,3}', str(seq)))
        print(f'Valid sequence: {seq}')

        return seq




def chkGC_content(sequence=Seq('GAGACCTT',generic_dna)):
    sequence = sequence.upper()
    seq_len = len(sequence)
    num_G, num_C = sequence.count('G'), sequence.count('C')
    GC_content = round((num_C + num_G) / seq_len * 100, 1)
    return GC_content




def GenOligoGC(length=25, GC_low=40, GC_high=60):
    generate = True
    #seed(int(time()))
    while generate == True:
        seq = randSeqGen(length)
        GC_content = GC(seq)
        GC_content1_10 = GC(seq[:10])
        GC_content11_20 = GC(seq[10:])
        if GC_content >= GC_low and GC_content <= GC_high:
            if GC_content1_10 >= GC_low and GC_content1_10 <= GC_high:
                if GC_content11_20 >= GC_low and GC_content11_20 <= GC_high:
                    return seq






def genPrimerPairs(primer_length=20, GC_low=40, GC_high=60):
    """Primer pairs for Half-Mers."""

    print('Primers for half-mers')

    forwTemplate5_3 = GenOligoGC(primer_length,GC_low, GC_high)
    """re.match checks if the first 2 Nuc are GC in the forward and backwards direction"""
    while not (re.match("[GC]{2}",str(forwTemplate5_3)) and
               re.match("[GC]{2}", str(forwTemplate5_3[::-1]))):

        forwTemplate5_3 = GenOligoGC(primer_length,GC_low, GC_high)

    forwPrimer =  forwTemplate5_3.reverse_complement()


    # print(f"Template   Seq 3\' - > 5\': {forwTemplate5_3[::-1]}")
    # print(f"ForwPrimer Seq 5\' - > 3\': {forwPrimer}")

    revPrimer = GenOligoGC(primer_length,GC_low, GC_high)
    """re.match checks if the first 2 Nuc are GC in the forward and backwards direction"""
    while not (re.match("[GC]{2}", str(revPrimer)) and
               re.match("[GC]{2}", str(revPrimer[::-1]))):

        revPrimer = GenOligoGC(primer_length, GC_low, GC_high)

    # print(f"RevPrimer  Seq 5\' - > 3\': {revPrimer}")

    return forwPrimer, revPrimer






def evalPrimerPairMT(fprimer, rprimer, ret_mt=False):
    """This will check the melting temperature

    The optimal melting temperature of the primers is 60–64°C, with
    an ideal temperature of 62°C, which is based on typical cycling and reaction conditions
    and the optimum temperature for PCR enzyme function. Ideally, the melting temperatures of
    the 2 primers should not differ by more than 2°C in order for both primers to bind
    simultaneously and efficiently amplify the product.
    PCR parameters used are from IDT: Oligo 0.2 uM Na 50 mM, Mg 3 mM, dNTPs 0.8 mM
    :param ret_mt: """

    fprimer_MT = round(MT.Tm_GC(fprimer, Na=50, Mg=3, dNTPs=0.8),2)
    rprimer_MT = round(MT.Tm_GC(rprimer, Na=50, Mg=3, dNTPs=0.8),2)

    fprimer_MT_NN = round(MT.Tm_NN(fprimer, Na=50, Mg=3, dNTPs=0.8, saltcorr=7),2)
    rprimer_MT_NN = round(MT.Tm_NN(fprimer, Na=50, Mg=3, dNTPs=0.8, saltcorr=7),2)


    print (f"forw primer: {fprimer}\nforw primer MT: {fprimer_MT} {fprimer_MT_NN} \n"
           f"rev  primer: {rprimer}\nrev primer MT : {rprimer_MT} {rprimer_MT_NN} \n")

    """Filters for primers that meet the MT standards"""
    if math.fabs(fprimer_MT - rprimer_MT) <= 3.0 and \
                max(fprimer_MT,rprimer_MT) < 65.0 and \
                min(fprimer_MT, rprimer_MT) > 59.0:

        print("MT of primer pair passed.\n")

        if ret_mt == False:
            return True
        else:
            return fprimer_MT, rprimer_MT
    else:
        print("MT for the primer pairs did not meet standards\n"); return False



def retPrimerPairMT(fprimer,rprimer):
    forwMT, revMT = [round(temp,2) for temp in evalPrimerPairMT(fprimer, rprimer, ret_mt=True)]
    return forwMT, revMT



def certMT_SelfDimerPrimerPairs(length=20, *, GC_low=40, GC_high=60, ret_str=True):

    chk_primer_mt = False
    while not chk_primer_mt:

        fprimer, rprimer = genPrimerPairs(primer_length=length, GC_low=GC_low, GC_high=GC_high)

        if selfDimerize(fprimer) and selfDimerize(rprimer):
            chk_primer_mt = evalPrimerPairMT(fprimer, rprimer)


    if ret_str == False:
        return fprimer,rprimer
    else:
        return str(fprimer), str(rprimer)





def selfDimerize(primer):
    if not isinstance(primer,Seq):
        primer = Seq(primer, generic_dna)
    length = len(primer)

    primer_Rev = primer[::-1]
    primer_Com = primer.complement()

    match = SequenceMatcher(a=primer_Rev, b=primer_Com).find_longest_match(0,length,0,length)

    all_match = SequenceMatcher(a=primer_Rev, b=primer_Com).get_matching_blocks()

    # print("rev  ", primer_Rev)
    # print('comp ', primer_Com)

    print(match)
    step = match.b - match.a
    print("step: ",step)
    print(all_match)
    phase_matches = []
    for match in all_match:
        diff = match.b - match.a
        if diff == step:
            phase_matches.append(match)
    phase_match_size = sum(match.size for match in phase_matches)
    # print(phase_matches)
    print("phase match: ",phase_match_size)

    if phase_match_size > 3:
        # print("rev  ",' ' * step, primer_Rev)
        # print('comp ',primer_Com)
        return False

    return True


def Orthogonality_Test(seq_1,seq_2,limit=4):
    if not isinstance(seq_1,Seq):
        seq_1 = Seq(seq_1, generic_dna)

    if not isinstance(seq_2,Seq):
        seq_2 = Seq(seq_2, generic_dna)


    if len(seq_1) == len(seq_2):
        length = len(seq_1)
    else:
        length = (max([len(seq_2),len(seq_2)]))

    print(f"Checking orthogonality between:\n{seq_1} and \n{seq_2}")

    seq_1_Rev = seq_1[::-1]
    seq_2_Com = seq_2.complement()

    match = SequenceMatcher(a=seq_1_Rev, b=seq_2_Com).find_longest_match(0, length, 0, length)

    all_match = SequenceMatcher(a=seq_1_Rev, b=seq_2_Com).get_matching_blocks()

    # print("rev  ", primer_Rev)
    # print('comp ', primer_Com)

    print(match)
    step = match.b - match.a
    print("step: ", step)
    print(all_match)
    phase_matches = []
    for match in all_match:
        diff = match.b - match.a
        if diff == step:
            phase_matches.append(match)
    phase_match_size = sum(match.size for match in phase_matches)
    # print(phase_matches)
    print("phase match: ", phase_match_size)

    if phase_match_size > limit:
        # print("rev  ",' ' * step, primer_Rev)
        # print('comp ',primer_Com)
        return False

    return True





def orthogonalPoolComb(tuple_list, len_output_seq=0, round=1):
    print(250 * "_")
    all_combinations = combinations(tuple_list,2)
    filtered_sequences = list()
    for pair1, pair2 in all_combinations:

        forwP_forwP = Orthogonality_Test(pair1.forw_primer,pair2.forw_primer,limit=5)
        forwP_revP = Orthogonality_Test(pair1.forw_primer,pair2.rev_primer,limit=5)
        revP_forwP = Orthogonality_Test(pair1.rev_primer, pair2.forw_primer,limit=5)
        revP_revP = Orthogonality_Test(pair1.rev_primer, pair2.rev_primer,limit=5)

        ortho_results = [forwP_forwP, forwP_revP, revP_forwP, revP_revP]

        if not (False in ortho_results):
            if pair1 not in filtered_sequences:
                filtered_sequences.append(pair1)
            if pair2 not in filtered_sequences:
                filtered_sequences.append(pair2)
        else:
            if pair1 in filtered_sequences:
                filtered_sequences.remove(pair1)
            if pair2 in filtered_sequences:
                filtered_sequences.remove(pair2)

    len_filtered_seq = len(filtered_sequences)
    if len_filtered_seq == len_output_seq and len_filtered_seq != 0:
        return filtered_sequences, round
    else:
        print(len_filtered_seq,round)
        return orthogonalPoolComb(filtered_sequences, len_output_seq=len_filtered_seq, round=round + 1)


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


        forwPrimer, revPrimer = certMT_SelfDimerPrimerPairs(length=length, GC_low=GC_low, GC_high=GC_high, ret_str=ret_str)
        if Orthogonality_Test(forwPrimer,revPrimer):
            primers[f'forw_primer_{ext}.{primerNum}'], primers[f'rev_primer_{ext}.{primerNum}'] = \
                forwPrimer,revPrimer
            primers[f'forw_MT_{ext}.{primerNum}'], primers[f'rev_MT_{ext}.{primerNum}'] = \
                retPrimerPairMT(primers[f'forw_primer_{ext}.{primerNum}'], primers[f'rev_primer_{ext}.{primerNum}'])

            pool[primerKeyName] = primers

        else:
            continue


    return pool





def gen5_3PrimerPairPools(pool_size=8,length=20,*,GC_low=40, GC_high=60,ret_str=True):
    """Return primer pair pools for the 5' extension and the 3' extension."""

    masterPool = dict()
    pool_5Ext = genPrimerPairPool(pool_size=pool_size,length=length,GC_low=GC_low,
                                  GC_high=GC_high,ret_str=ret_str)

    pool_3Ext = genPrimerPairPool(pool_size=pool_size, length=length, GC_low=GC_low,
                                  GC_high=GC_high, ret_str=ret_str, ext=3)

    masterPool["5\'_ext_pool"] = pool_5Ext; masterPool["3\'_ext_pool"] = pool_3Ext

    return masterPool






def retFormattedSeq(seq,rand_reg_len=15):
    ret_list =list()
    seqTemplate = namedtuple('SeqTemplate',['pair_num','extension','forw_primer','forw_MT', 'rev_primer','rev_MT',
                                             'template','sequence'])
    for primer_pair in seq:
        forw_primer = Seq(primer_pair.forw_primer,alphabet=generic_dna)
        forwP_template = forw_primer.reverse_complement()
        seq_template = primer_pair.rev_primer + '-'*rand_reg_len + forwP_template
        ret_list.append(seqTemplate(*primer_pair,str(forwP_template),str(seq_template)))
    # pprint(ret_list)
    # print('\n\n')
    return ret_list







def retLigatedSeqData(*,filename,five_prime,three_prime, ret_dict=True):
    five_prime= 'Pair_'+ str(five_prime); three_prime = "Pair_" + str(three_prime)
    SeqData = namedtuple('LigSeqData', ['pair_num_5', 'pair_num_3', 'template','forw_primer', 'rev_primer',
                                             'sequence','oligo_bridge'])

    with open(filename,'r') as f:
        file = json.load(f)
        five_prime_data = file.get(five_prime)
        three_prime_data = file.get(three_prime)
        five_prime_seq = five_prime_data['sequence']
        three_prime_seq = three_prime_data['sequence']


        #return five_prime_data, three_prime_data

        print('five_prime_extension', five_prime_data['sequence'], five_prime_data['sequence'][:5],sep='\t\t')
        print('three_prime_extension', three_prime_data['sequence'], three_prime_data['sequence'][-5:],sep='\t\t')

    five_prime_extension = five_prime_seq[:10]
    three_prime_extension = three_prime_seq[-10:]
    ligated_extensions = three_prime_extension + five_prime_extension
    ligated_sequence = three_prime_seq + "-" + five_prime_seq
    ligated_extensions = Seq(ligated_extensions)
    oligo_bridge= ligated_extensions.reverse_complement()[5:15]
    oligo_bridge_var = oligo_bridge[:5] + "-" + oligo_bridge[5:]

    final_data = SeqData(five_prime_data['pair_num'], three_prime_data['pair_num'], five_prime_data['template'],
                   five_prime_data['forw_primer'], three_prime_data['rev_primer'],
                   ligated_sequence, str(oligo_bridge_var))
    if ret_dict:
        return final_data._asdict()
    else:
        return final_data






def NameTupleAllPrimers(pool_dict=None):
    if pool_dict == None:
        pool_dict = gen5_3PrimerPairPools(10)

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





def export2File(dict_, *, filename, f_ext='json'):
    now = datetime.now().strftime("%m-%d-%y_%H:%M")
    with open(f"{filename}-{now}.{f_ext}", 'w+') as f:
        json.dump(dict_, f, indent=4, sort_keys=True)
    print('created file')





def convNamedTuple2Dict(namedtuples):
    ret_dict = dict()
    for index, ntuple in enumerate(namedtuples):
        print(type(ntuple))
        print(ntuple)
        ret_dict[f'Pair_{index + 1}'] = ntuple._asdict()

    return ret_dict




def grabPrimersFromFile(filename, *, chosen_3_pairs=None, chosen_5_pairs=None):
    if chosen_3_pairs == None:
        chosen_3_pairs = ['Pair_21']
    if chosen_5_pairs == None:
        chosen_5_pairs = ['Pair_24']


    chosen_primers = list()
    with open(filename, 'r') as f:
        file = json.load(f)
        PrimerPairs = namedtuple('PrimerPair',['pair_num','extension','forw_primer','forw_MT', 'rev_primer','rev_MT'])
        for index,(three_prime, five_prime) in enumerate(zip(chosen_3_pairs,chosen_5_pairs)):
            chosen_primers.append(PrimerPairs(three_prime,**file[three_prime]))
            chosen_primers.append(PrimerPairs(five_prime,**file[five_prime]))

    return chosen_primers


def genSplintSeq(filename, *, splint_len=20):
    chosen_3_pairs = ['Pair_2', 'Pair_12', 'Pair_36', 'Pair_39', 'Pair_45']
    chosen_5_pairs = ['Pair_14', 'Pair_17', 'Pair_19', 'Pair_29', 'Pair_35']
    chosen_primers = list()
    with open(filename, 'r') as f:
        file = json.load(f)

    half_splint_len = splint_len//2
    splint_oligo_complement = ''
    for dic in file.values():
        if dic['extension'] == '3_prime':
            oligo_4_splint = dic['sequence'][-half_splint_len:]
            splint_oligo_complement = oligo_4_splint + splint_oligo_complement

        elif dic['extension'] == '5_prime':
            oligo_4_splint = dic['sequence'][:half_splint_len]
            splint_oligo_complement += oligo_4_splint
    splint_oligo_complement = Seq(splint_oligo_complement)
    print("Sequence:",splint_oligo_complement,
          "\nReverse Compliment:", splint_oligo_complement.reverse_complement(),end='\n\n')

    return splint_oligo_complement.reverse_complement()








if __name__=="__main__":

    #print([certMTPrimerPairs() for _ in range(5)])
    #print(certMTPrimerPairs())

    # pool = gen5_3PrimerPairPools(20)
    # print(pool)
    # pool_tuple = NameTupleAllPrimers(pool)
    # print(pool_tuple)
    # print(250 * "_")
    # ortho_tuple, round = orthogonalPoolComb(pool_tuple)
    # print(ortho_tuple)
    # pool_dict = convNamedTuple2Dict(ortho_tuple)
    # export2File(pool_dict,filename='OrthoPrimerPairs',f_ext='txt')

    # chosen_primers = grabPrimersFromFile("OrthoPrimerPairs-11-18-18_22:45.txt")
    # print(chosen_primers)
    # formatted_seq = retFormattedSeq(chosen_primers)
    # formatted_dict = convNamedTuple2Dict(formatted_seq)
    # export2File(formatted_dict, filename='ChosenPrimerPairs', f_ext='txt')


    # final_seq = retLigatedSeqData(filename="ChosenPrimerPairs-11-19-18_06:00.txt",five_prime=21,three_prime=24)
    # print(final_seq)


    genSplintSeq('ChosenPrimerPairs-11-19-18_06:00.txt')
    genSplintSeq('ChosenPrimerPairs-11-19-18_06:00.txt',splint_len=16)
    genSplintSeq('ChosenPrimerPairs-11-19-18_06:00.txt', splint_len=12)

