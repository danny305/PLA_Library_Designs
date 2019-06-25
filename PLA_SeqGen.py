from numpy.random import choice
from itertools import combinations, product
from difflib import SequenceMatcher

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import GC, MeltingTemp as MT

from math import fabs
from datetime import datetime

import re
import json


class PLA_Seq():

    def __init__(self,orth_seqs=[], seq_len=20,):

        if not isinstance(orth_seqs, (list,tuple)):
            raise TypeError('Parameter orth_seq must be of type List or Tuple.')
        self.ref_orth_sequence =tuple(orth_seqs) if len(orth_seqs) > 0 else None
        self.seq_len = seq_len


    @staticmethod
    def randSeqGen(length=20, clamps=True):
        """
        Generates and returns a DNA Sequence object of the desired nucleotide length.
        Where the maximum number of adjacent same nucleotides is 3
        """

        seq_approved = False

        p_length = length - 6 if clamps else length

        def return_rand_clamp(end):
            """
            :param end: int 5 or 3 should be passed. This specifies whether its the 5 prime or
                3 prime clamp.
            :return: 3 nucleotide clamp.
            """
            if end is 5:
                combos = tuple(''.join(pair)
                               for pair in combinations(['CC', 'CG', 'GC', 'GG', 'A', 'T'], 2)
                               if len(pair[0] + pair[1]) == 3)
            elif end is 3:
                combos = tuple(''.join(pair)
                               for pair in combinations(['A', 'T', 'CC', 'CG', 'GC', 'GG'], 2)
                               if len(pair[0] + pair[1]) == 3)
            else:
                raise AssertionError('A 5 or 3 are the designed input parameter.')

            return choice(combos)


        seq = return_rand_clamp(5) if clamps else ""

        while not seq_approved:

            test_seq = seq + ''.join([choice(['A', 'T', 'C', 'G']) for _ in range(p_length)])

            if clamps:
                test_seq += return_rand_clamp(3)

            if re.search("[GC]{4,100}|[G]{4,100}|[C]{4,100}|[A]{4,100}|[T]{4,100}", test_seq):
                #print(f"The sequence: {test_seq} is invalid")
                continue

            if len(re.findall('C{2,3}|G{2,3}|A{2,3}|T{2,3}', test_seq)) > 3 or \
                    len(re.findall('C{3}|G{3}|A{3}|T{3}', test_seq)) > 1:
                #print(f'sequence not random enough: {test_seq}')
                continue


            seq_approved = True
            seq = test_seq
            #print(re.findall('[C]{2,3}|[G]{2,3}|[A]{2,3}|[T]{2,3}', seq))

        #print(f'Valid sequence: {seq}')
        return Seq(seq, generic_dna)



    @staticmethod
    def chkGC_content(sequence=Seq('GAGACCTT', generic_dna)):
        sequence = sequence.upper()
        seq_len = len(sequence)
        num_G, num_C = sequence.count('G'), sequence.count('C')
        GC_content = round((num_C + num_G) / seq_len * 100, 1)
        return GC_content



    @staticmethod
    def GenOligoGC(length=20, GC_low=40, GC_high=60, clamps=True):
        while True:
            seq = PLA_Seq.randSeqGen(length, clamps=clamps)
            GC_content = GC(seq)
            GC_content_f_half = GC(seq[:length//2])
            GC_content_b_half = GC(seq[length//2:])
            if GC_content >= GC_low and GC_content <= GC_high:
                if GC_content_f_half >= GC_low and GC_content_f_half <= GC_high:
                    if GC_content_b_half >= GC_low and GC_content_b_half <= GC_high:
                        #print(f'Valid GC content: {seq}  GC%: {GC_content}',)
                        return seq



    @staticmethod
    def genPrimerPairs(primer_length=20, GC_low=40, GC_high=60, clamps=True):
        """Primer pairs for Half-Mers."""
        #print('Primers for half-mers')

        forwTemplate5_3 = PLA_Seq.GenOligoGC(primer_length, GC_low, GC_high, clamps)
        forwPrimer = forwTemplate5_3.reverse_complement()
        revPrimer = PLA_Seq.GenOligoGC(primer_length, GC_low, GC_high, clamps)

        # print(f"Template   Seq 3\' - > 5\': {forwTemplate5_3[::-1]}")
        # print(f"ForwPrimer Seq 5\' - > 3\': {forwPrimer}")
        # print(f"RevPrimer  Seq 5\' - > 3\': {revPrimer}")

        return forwPrimer, revPrimer



    @staticmethod
    def evalPrimerPairMT(fprimer, rprimer, ret_mt=False):
        """This will check the melting temperature

        The optimal melting temperature of the primers is 60–64°C, with
        an ideal temperature of 62°C, which is based on typical cycling and reaction conditions
        and the optimum temperature for PCR enzyme function. Ideally, the melting temperatures of
        the 2 primers should not differ by more than 2°C in order for both primers to bind
        simultaneously and efficiently amplify the product.
        PCR parameters used are from IDT: Oligo 0.2 uM Na 50 mM, Mg 3 mM, dNTPs 0.8 mM
        :param ret_mt: """

        fprimer_MT = round(MT.Tm_GC(fprimer, Na=50, Mg=3, dNTPs=0.8), 2)
        rprimer_MT = round(MT.Tm_GC(rprimer, Na=50, Mg=3, dNTPs=0.8), 2)

        fprimer_MT_NN = round(MT.Tm_NN(fprimer, Na=50, Mg=3, dNTPs=0.8, saltcorr=7), 2)
        rprimer_MT_NN = round(MT.Tm_NN(fprimer, Na=50, Mg=3, dNTPs=0.8, saltcorr=7), 2)

        #print(f"forw primer: {fprimer}\nforw primer MT: {fprimer_MT} (NN){fprimer_MT_NN} \n"
         #     f"rev  primer: {rprimer}\nrev primer MT : {rprimer_MT} (NN){rprimer_MT_NN} \n")

        """Filters for primers that meet the MT standards"""
        if fabs(fprimer_MT - rprimer_MT) <= 3.0 and \
                max(fprimer_MT, rprimer_MT) < 65.0 and \
                min(fprimer_MT, rprimer_MT) > 59.0:

            #print("MT of primer pair passed.\n")

            if ret_mt == False:
                return True
            else:
                return fprimer_MT, rprimer_MT
        else:
            #print("MT for the primer pairs did not meet standards\n")
            return False



    @staticmethod
    def evalSeqMT(seq, ret_mt=False):
        seq_MT = round(MT.Tm_GC(seq, Na=50, Mg=3, dNTPs=0.8), 2)

        if seq_MT < 65.0 and seq_MT > 59.0:
            if ret_mt == False:
                return True
            else:
                return seq_MT
        else:
            # print("MT for the primer pairs did not meet standards\n")
            return False




    @staticmethod
    def retPrimerPairMT(fprimer, rprimer):
        forwMT, revMT = [round(temp, 2)
                            for temp in PLA_Seq.evalPrimerPairMT(fprimer, rprimer, ret_mt=True)]
        return forwMT, revMT




    @staticmethod
    def certMT_SelfDimerPrimerPairs(length=20, *, GC_low=40, GC_high=60, clamps=True, ret_str=True):

        chk_primer_mt = False
        while not chk_primer_mt:
            fprimer, rprimer = PLA_Seq.genPrimerPairs(primer_length=length,
                                                      GC_low=GC_low, GC_high=GC_high, clamps=clamps)

            if PLA_Seq.selfDimerizeTest(fprimer) and PLA_Seq.selfDimerizeTest(rprimer):
                chk_primer_mt = PLA_Seq.evalPrimerPairMT(fprimer, rprimer)

        return (fprimer, rprimer) if ret_str == False else (str(fprimer), str(rprimer))




    @staticmethod
    def genCertMT_SelfDimerSeq(length=20, *, GC_low=40, GC_high=60, clamps=True, ret_str=True):

        chk_primer_mt = False
        while not chk_primer_mt:
            seq1 = PLA_Seq.GenOligoGC(length, GC_low, GC_high, clamps)

            if PLA_Seq.selfDimerizeTest(seq1):
                chk_primer_mt = PLA_Seq.evalSeqMT(seq1)

        return seq1 if ret_str == False else str(seq1)





    @staticmethod
    def calcPhaseMatch(seq1, seq2, limit = 3):

        if not (isinstance(seq1, Seq) and isinstance(seq2, Seq)):
            raise AssertionError('Seq 1 and Seq 2 must both be either a str or Seq object.')

        # if len(seq1) != len(seq2):
        #     raise AssertionError('Seq 1 and Seq 2 must be of the same length.')

        length = len(seq1) if len(seq1) == len(seq2) else max([len(seq1), len(seq2)])
        match_obj = SequenceMatcher(a=seq1, b=seq2)
        longest_match = match_obj.find_longest_match(0, length, 0, length)
        if longest_match.size >3:
            return False
        all_match = match_obj.get_matching_blocks()
        phase_gap = longest_match.b - longest_match.a

        # print("seq1  ", seq1)
        # print('seq2  ', seq2)
        # print('Match: ', longest_match)
        # print("phase_gap: ", phase_gap)
        # print("all_match", all_match)

        phase_match_size = sum([match.size
                                    for match in all_match
                                        if (match.b - match.a) == phase_gap])

        # phase_matches = []
        # for index, match in enumerate(all_match):
        #     diff = match.b - match.a
        #     if diff == phase_gap:
        #         phase_matches.append(match)
        #         print("index: ", index)
        #         size = match.size
        #         print("rev  ", seq1[match.a:match.a + size] + " " + seq1[match.a + size:])
        #         print('comp ', seq2[match.b:match.b + size] + " " + seq2[match.b + size:], end='\n\n')
        # phase_match_size = sum(match.size for match in phase_matches)
        # print(phase_matches)

        #print("phase match size: ", phase_match_size)
        #print("phase match size: ", phase_match_size_1)

        # if phase_match_size <=limit:
        #     print(phase_match_size)
        #     return True
        # return False

        return True if phase_match_size <= limit else False






    @staticmethod
    def selfDimerizeTest(primer):
        if not isinstance(primer, Seq):
            primer = Seq(primer, generic_dna)
        length = len(primer)

        primer_Rev = primer[::-1]
        primer_Com = primer.complement()

        return PLA_Seq.calcPhaseMatch(primer_Rev, primer_Com)





    @staticmethod
    def OrthogonalityTest(seq_1, seq_2, limit=3):
        if isinstance(seq_1, str) and isinstance(seq_2, str):
            seq_1 = Seq(seq_1, generic_dna)
            seq_2 = Seq(seq_2, generic_dna)

        elif not(isinstance(seq_1, Seq) and isinstance(seq_2, Seq)):
            raise AssertionError('Seq 1 and Seq 2 must both be either a str or Seq object.')

        #print(f"Checking orthogonality between:\n{seq_1} and \n{seq_2}")
        seq_1_Rev = seq_1[::-1]
        seq_2_Com = seq_2.complement()

        return PLA_Seq.calcPhaseMatch(seq_1_Rev, seq_2_Com, limit)



    @staticmethod
    def genOrthPrimerPairs(length=20, limit=4, GC_low=40, GC_high=60, ret_str=True, clamps=True):

        gen_orthogonal_pair = False
        while not gen_orthogonal_pair:
            forwPrimer, revPrimer = PLA_Seq.certMT_SelfDimerPrimerPairs(
                length=length, GC_low=GC_low, GC_high=GC_high, ret_str=ret_str, clamps=clamps)

            if PLA_Seq.OrthogonalityTest(forwPrimer, revPrimer, limit):
                return forwPrimer, revPrimer



    @staticmethod
    def genOrthoPrimerPairPool(pool_size=4, length=20, *, GC_low=40, GC_high=60, ret_str=True, limit=4, clamps=True):
        """Generate a pool of primer pairs that pass the melt temp criteria defined in
        genCertPrimerPairs(). Default is the 5'extension half-asstemers primers."""

        pool = dict()
        count = 1

        pool['meta_data'] = {
            'pool_size': pool_size, 'primer_length': length,
            "GC_low": GC_low, "GC_high": GC_high,
            "date": datetime.today().strftime("%-m-%d-%Y_%-H:%-M")
        }

        while len(pool.keys()) <= pool_size:
            primerKeyName = f"Pair-{count}"
            primers = dict()

            forwPrimer, revPrimer = PLA_Seq.genOrthPrimerPairs(length=length, GC_low=GC_low, GC_high=GC_high,
                                                               ret_str=ret_str,limit=limit, clamps=clamps)

            if PLA_Seq.poolOrthogonalityTest(pool, forwPrimer, revPrimer, limit):
                primers['Seq_1'], primers['Seq_2'] = forwPrimer, revPrimer
                primers['Seq_1_MT'], primers['Seq_2_MT'] = PLA_Seq.retPrimerPairMT(forwPrimer, revPrimer)
                pool[primerKeyName] = primers
                print(f'Added primer pair {count}')
                count += 1

        return pool




    def genOrthoPrimers2Ref(self, pool_size=3, length=20, *, GC_low=40, GC_high=60, ret_str=True,
                            limit=4, filename='OrthoPool_OrthoRef', clamps=True, chk_orth_2_new_seq=True):

        if len(self.ref_orth_sequence) == 0:
            return PLA_Seq.genOrthoPrimerPairPool(pool_size=pool_size, length=length,
                                                  GC_low=GC_low, GC_high=GC_high,
                                                  ret_str=ret_str, limit=limit)

        pool = dict()
        count = 1
        attempt = 0

        pool['meta_data'] = {
            'pool_size': pool_size, 'primer_length': length,
            "GC_low": GC_low, "GC_high": GC_high,
            "date": datetime.today().strftime("%-m-%d-%Y_%-H:%-M"),
            "ortho_limit": limit,
            'reference_sequences': self.ref_orth_sequence,
            'seq_ortho_2_each_other': chk_orth_2_new_seq,
            'seq_ortho_2_ref': True
        }

        while len(pool.keys()) <= pool_size:
            attempt += 1
            primerKeyName = f"Pair-{count}"
            primers = dict()

            forwPrimer, revPrimer = PLA_Seq.genOrthPrimerPairs(length=length, GC_low=GC_low, GC_high=GC_high,
                                                               ret_str=ret_str,limit=limit, clamps=clamps)

            # combos = [tuple(combinations([forwPrimer, revPrimer, ref_seq],2))[1],
            #             for ref_seq in self.ref_orth_sequence]

            combos = tuple(product([forwPrimer, revPrimer], self.ref_orth_sequence, repeat=1))
            if not all([PLA_Seq.OrthogonalityTest(*ppair, limit=limit) for ppair in combos]):
                print(f'{attempt} Failed orthogonality to current primers - forw: {forwPrimer}   rev: {revPrimer}')
                continue
            print(f'Forw: {forwPrimer} Rev: {revPrimer} are orthogonal to current primers.')



            if chk_orth_2_new_seq and PLA_Seq.poolOrthogonalityTest(pool, forwPrimer, revPrimer, limit):
                primers['Seq_1'], primers['Seq_2'] = forwPrimer, revPrimer
                primers['Seq_1_MT'], primers['Seq_2_MT'] = PLA_Seq.retPrimerPairMT(forwPrimer, revPrimer)
                pool[primerKeyName] = primers
                print(f'Added primer pair {count}')
                self.ortho_pool_ortho_ref = pool
                self.exportPool2_JSON(filename)
                count += 1

            else:
                primers['Seq_1'], primers['Seq_2'] = forwPrimer, revPrimer
                primers['Seq_1_MT'], primers['Seq_2_MT'] = PLA_Seq.retPrimerPairMT(forwPrimer, revPrimer)
                pool[primerKeyName] = primers
                print(f'Added primer pair {count}')
                self.ortho_pool_ortho_ref = pool
                self.exportPool2_JSON(filename)
                count += 1



        self.ortho_pool_ortho_ref = pool
        return pool



    def genOrthoSeqPool2Ref(self, pool_size=3, length=20, *, GC_low=40, GC_high=60, ret_str=True,
                            limit=6, filename='OrthoPool_OrthoRef', clamps=True, chk_orth_2_new_seq=True):

        if len(self.ref_orth_sequence) == 0:
            raise ValueError("Missing reference sequences to compare new sequences to.")

        pool = dict()
        count = 1
        attempt = 0

        pool['meta_data'] = {
            'pool_size': pool_size, 'primer_length': length,
            "GC_low": GC_low, "GC_high": GC_high,
            "date": datetime.today().strftime("%-m-%d-%Y_%-H:%-M"),
            "ortho_limit": limit,
            'reference_sequences': self.ref_orth_sequence,
            'seq_ortho_2_each_other': chk_orth_2_new_seq,
            'seq_ortho_2_ref': True
        }

        while len(pool.keys()) <= pool_size:
            attempt += 1
            keyName = f"Ortho-Seq-{count}"
            sequence = dict()

            seq1 = PLA_Seq.genCertMT_SelfDimerSeq(length=length, GC_low=GC_low, GC_high=GC_high,
                                                               ret_str=ret_str, clamps=clamps)

            seq2 = PLA_Seq.genCertMT_SelfDimerSeq(length=length, GC_low=GC_low, GC_high=GC_high,
                                                  ret_str=ret_str, clamps=clamps)

            # combos = [tuple(combinations([forwPrimer, revPrimer, ref_seq],2))[1],
            #             for ref_seq in self.ref_orth_sequence]

            combos1 = tuple(product([seq1], self.ref_orth_sequence, repeat=1))
            combos2 = tuple(product([seq2], self.ref_orth_sequence, repeat=1))


            if all([PLA_Seq.OrthogonalityTest(*ppair, limit=limit) for ppair in combos1]):
                if chk_orth_2_new_seq and PLA_Seq.poolOrthogonalityTest(pool, seq1, limit):
                    sequence['Seq_1'] = seq1
                    sequence['Seq_1_MT'] = PLA_Seq.evalSeqMT(seq1, ret_mt=True)
                    pool[keyName] = sequence
                    print(f'Seq1: {seq1} is orthogonal to current sequence.')
                    self.ortho_pool_ortho_ref = pool
                    self.exportPool2_JSON(filename)
                    print(f'Added sequence: {count}')
                    count += 1


                else:
                    sequence['Seq_1']  = seq1
                    sequence['Seq_1_MT'] = PLA_Seq.evalSeqMT(seq1, ret_mt=True)
                    pool[keyName] = sequence
                    print(f'Seq1: {seq1} is orthogonal to current sequence.')
                    self.ortho_pool_ortho_ref = pool
                    self.exportPool2_JSON(filename)
                    print(f'Added sequence: {count}')
                    count += 1
            else:
                pass
                #print('Combos 1: ', [PLA_Seq.OrthogonalityTest(*ppair, limit=limit) for ppair in combos1])




            if all([PLA_Seq.OrthogonalityTest(*ppair, limit=limit) for ppair in combos2]):
                if chk_orth_2_new_seq and PLA_Seq.poolOrthogonalityTest(pool, seq2, limit):
                    sequence['Seq_1'] = seq2
                    _, sequence['Seq_2_MT'] = PLA_Seq.evalSeqMT(seq2, ret_mt=True)
                    pool[keyName] = sequence
                    print(f'Seq2: {seq2} is orthogonal to current sequence.')
                    self.ortho_pool_ortho_ref = pool
                    self.exportPool2_JSON(filename)
                    print(f'Added sequence: {count}')
                    count += 1
                    continue

                else:
                    sequence['Seq_2']  = seq2
                    sequence['Seq_2_MT'] = PLA_Seq.evalSeqMT(seq2, ret_mt=True)
                    pool[keyName] = sequence
                    print(f'Seq2: {seq2} is orthogonal to current sequence.')
                    self.ortho_pool_ortho_ref = pool
                    self.exportPool2_JSON(filename)
                    print(f'Added sequence: {count}')
                    count += 1
                    continue
            else:
                pass
                #print('Combos 2: ', [PLA_Seq.OrthogonalityTest(*ppair, limit=limit) for ppair in combos1])


            print(f'{attempt} Failed orthogonality to current sequences - Seq1: {seq1}   Seq2: {seq2}')

        self.ortho_pool_ortho_ref = pool
        return pool



    def exportPool2_JSON(self, filename='OrthoPool_OrthoRef'):
        #now = datetime.now().strftime("%m-%d-%y_%H:%M")
        now = self.ortho_pool_ortho_ref['meta_data']['date']
        name = f"./pla_primer_output/{filename}-{now}.json"
        with open(name, 'w+') as f:
            json.dump(self.ortho_pool_ortho_ref, f, indent=4, sort_keys=True)
        print(f'created file: {name}')


    @staticmethod
    def poolOrthogonalityTest(pool, forwPrimer, revPrimer, limit=3):

        for key, p_dict in pool.items():
            if key is 'meta_data':
                continue

            # print([PLA_Seq.OrthogonalityTest(*ppair, limit=4) for ppair in combos])
            # print(all([PLA_Seq.OrthogonalityTest(*ppair, limit=4) for ppair in combos]))

            combos = list(combinations([p_dict['Seq_1'], p_dict['Seq_2'], forwPrimer], 2))[1:] + \
                     list(combinations([p_dict['Seq_1'], p_dict['Seq_2'], revPrimer], 2))[1:]

            if not all([PLA_Seq.OrthogonalityTest(*ppair, limit=limit) for ppair in combos]):
                return False
        return True



    @staticmethod
    def export2File(dict_, *, filename, f_ext='json'):
        now = datetime.now().strftime("%m-%d-%y_%H:%M")
        with open(f"{filename}-{now}.{f_ext}", 'w+') as f:
            json.dump(dict_, f, indent=4, sort_keys=True)
        print('created file')



    @staticmethod
    def genSplintSeq(seq_3E, seq_5E, length=[12,16,20]):


        return [(splint_len, PLA_Seq.retComplementStrand(seq_3E[-1*(splint_len//2):] + seq_5E[:splint_len//2]))
                    for splint_len in length]


        # print(f'3E sequence: {seq_3E}')
        # print(f'5E sequence: {seq_5E}')
        # for splint_len in length:
        #     splint_template = seq_3E[-1*(splint_len//2):] + seq_5E[:splint_len//2]
        #     splint_seq = PLA_Seq.retComplementStrand(splint_template)
        #
        #     print(f"Splint Template: {splint_template} length: {len(splint_template)}" )
        #     print(f'Splint Sequence: {splint_seq} length: {len(splint_seq)}')



    @staticmethod
    def retComplementStrand(strand):
        if isinstance(strand, str):
            strand = Seq(strand, generic_dna)
        elif not isinstance(strand, Seq):
            raise TypeError('strand argument must be either a str or Seq type.')

        return strand.reverse_complement()







def main():
    # primer_pairs = PLA_Seq.genPrimerPairs()
    # PLA_Seq.selfDimerizeTest(primer_pairs[0])
    # PLA_Seq.evalPrimerPairMT(*primer_pairs)
    #print(PLA_Seq.certMT_SelfDimerPrimerPairs())

    #print(PLA_Seq.genOrthoPrimerPairPool(3))
    #print(PLA_Seq.genOrthPrimerPairs())

    # Generate orthogonal sequence to current primers being used in lab
    # obj = PLA_Seq(['GGTGTAAAGTCCACTCTACC', 'GCTCAGACCAATGGAGATGC',
    #                'CCAGTCAGTAGTAACGCTGC', 'GGTCAGTGTTGGATACGAGC'])
    # obj.genOrthoPrimers2Ref()


    # Generate orthogonal sequence for middle regions to new primers being purchased
    # obj = PLA_Seq(['CGAGATGGAGAAACAGGAGC', 'CGAAGTAGGGATTAGTGTGC',
    #                'CCAGTCTACGAGTCTTGTCC', 'GCTAACAGGCACGCAGTTGC'])
    # obj.genOrthoPrimers2Ref(filename='NewPrimers_ortho_IntraSeq', clamps=False)



    # Generate orthogonal 3E forw primer
    # obj = PLA_Seq([#'GGTGTAAAGTCCACTCTACC', 'GGTCAGTGTTGGATACGAGC',
    #                'CCAGTCAGTAGTAACGCTGC',  'GCTCAGACCAATGGAGATGC',
    #                'CGAAGTAGGGATTAGTGTGC', 'TGTCAGTGGTTCAGTCAGCA',
    #                'GGACAAGACTCGTAGACTGG', 'CCAGTCTACGAGTCTTGTCC',
    #                'GCTAACAGGCACGCAGTTGC', 'AGTAAGTGAGCAGTGGGTTG'])
    # obj.genOrthoSeqPool2Ref(pool_size=8,filename='3E_ortho_forwPrimer', clamps=True,
    #                         chk_orth_2_new_seq=False, limit=5)



    # obj = PLA_Seq([  # 'GGTGTAAAGTCCACTCTACC', 'GGTCAGTGTTGGATACGAGC',
    #     'CCAGTCAGTAGTAACGCTGC', 'GCTCAGACCAATGGAGATGC',])
    # obj.genOrthoSeqPool2Ref(filename='3E_ortho_forwPrimer', clamps=True,
    #                     chk_orth_2_new_seq=False, limit=5)

    # Generate orthogonal 5E rev primer
    # obj = PLA_Seq(['GGTGTAAAGTCCACTCTACC', 'GCTCAGACCAATGGAGATGC',
    #                'CCAGTCAGTAGTAACGCTGC', 'GGTCAGTGTTGGATACGAGC',
    #                'GCTCCTGTTTCTCCATCTCG', 'CGAGATGGAGAAACAGGAGC',
    #                'CGAAGTAGGGATTAGTGTGC', 'TGTCAGTGGTTCAGTCAGCA',
    #                'GGACAAGACTCGTAGACTGG', 'CCAGTCTACGAGTCTTGTCC',
    #                'AGTAAGTGAGCAGTGGGTTG'
    #                ])
    # obj.genOrthoPrimers2Ref(filename='3E_ortho_forwPrimer', clamps=True,
    #                         chk_orth_2_new_seq = False)



    print(PLA_Seq.genSplintSeq('GCTTTCTCTGATACCTGACG','GCTAACAGGCACGCAGTTGC'))









# obj1 = PLA_Seq()
# obj2 = PLA_Seq(orth_seqs=['Hello'])
# obj3 = PLA_Seq()
#
# obj4 = PLA_Seq(orth_seqs=('Hello',))

#print(PLA_Seq.genPrimerPairs())

if __name__=="__main__":
    main()


