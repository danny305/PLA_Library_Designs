from numpy.random import choice, seed
from itertools import combinations
from difflib import SequenceMatcher

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import GC, MeltingTemp as MT

from math import fabs

import re


class PLA_Seq():

    def __init__(self,orth_seqs=[], seq_len=20,):

        if not isinstance(orth_seqs, (list,tuple)):
            raise TypeError('Parameter orth_seq must be of type List or Tuple.')
        self.input_orth_sequence =tuple(orth_seqs) if len(orth_seqs) > 0 else None
        self.seq_len = seq_len


    @staticmethod
    def randSeqGen(length=20):
        """
        Generates and returns a DNA Sequence object of the desired nucleotide length.
        Where the maximum number of adjacent same nucleotides is 3
        """

        seq_approved = False

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

        seq = return_rand_clamp(5)

        while not seq_approved:
            test_seq = seq + ''.join([choice(['A', 'T', 'C', 'G']) for _ in range(length - 6)])


            if re.search("[G]{4,100}|[C]{4,100}|[A]{4,100}|[T]{4,100}", test_seq):
                #print(f"The sequence: {test_seq} is invalid")
                continue

            if len(re.findall('C{2,3}|G{2,3}|A{2,3}|T{2,3}', test_seq)) > 3 or \
                    len(re.findall('C{3}|G{3}|A{3}|T{3}', test_seq)) > 1:
                #print(f'sequence not random enough: {test_seq}')
                continue


            seq_approved = True
            seq = test_seq
            #print(re.findall('[C]{2,3}|[G]{2,3}|[A]{2,3}|[T]{2,3}', seq))


        seq += return_rand_clamp(3)
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
    def GenOligoGC(length=20, GC_low=40, GC_high=60):
        generate = True
        # seed(int(time()))
        while generate == True:
            seq = PLA_Seq.randSeqGen(length)
            GC_content = GC(seq)
            GC_content_f_half = GC(seq[:length//2])
            GC_content_b_half = GC(seq[length//2:])
            if GC_content >= GC_low and GC_content <= GC_high:
                if GC_content_f_half >= GC_low and GC_content_f_half <= GC_high:
                    if GC_content_b_half >= GC_low and GC_content_b_half <= GC_high:
                        print(f'Valid GC content: {seq}  GC%: {GC_content}',)
                        return seq


    @staticmethod
    def genPrimerPairs(primer_length=20, GC_low=40, GC_high=60):
        """Primer pairs for Half-Mers."""
        #print('Primers for half-mers')

        forwTemplate5_3 = PLA_Seq.GenOligoGC(primer_length, GC_low, GC_high)
        forwPrimer = forwTemplate5_3.reverse_complement()
        revPrimer = PLA_Seq.GenOligoGC(primer_length, GC_low, GC_high)

        print(f"Template   Seq 3\' - > 5\': {forwTemplate5_3[::-1]}")
        print(f"ForwPrimer Seq 5\' - > 3\': {forwPrimer}")
        print(f"RevPrimer  Seq 5\' - > 3\': {revPrimer}")

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

        print(f"forw primer: {fprimer}\nforw primer MT: {fprimer_MT} (NN){fprimer_MT_NN} \n"
              f"rev  primer: {rprimer}\nrev primer MT : {rprimer_MT} (NN){rprimer_MT_NN} \n")

        """Filters for primers that meet the MT standards"""
        if fabs(fprimer_MT - rprimer_MT) <= 3.0 and \
                max(fprimer_MT, rprimer_MT) < 65.0 and \
                min(fprimer_MT, rprimer_MT) > 59.0:

            print("MT of primer pair passed.\n")

            if ret_mt == False:
                return True
            else:
                return fprimer_MT, rprimer_MT
        else:
            print("MT for the primer pairs did not meet standards\n");
            return False


    @staticmethod
    def retPrimerPairMT(fprimer, rprimer):
        forwMT, revMT = [round(temp, 2)
                            for temp in PLA_Seq.evalPrimerPairMT(fprimer, rprimer, ret_mt=True)]
        return forwMT, revMT



    @staticmethod
    def certMT_SelfDimerPrimerPairs(length=20, *, GC_low=40, GC_high=60, ret_str=True):

        chk_primer_mt = False
        while not chk_primer_mt:
            fprimer, rprimer = PLA_Seq.genPrimerPairs(primer_length=length, GC_low=GC_low, GC_high=GC_high)

            if PLA_Seq.selfDimerizeTest(fprimer) and PLA_Seq.selfDimerizeTest(rprimer):
                chk_primer_mt = PLA_Seq.evalPrimerPairMT(fprimer, rprimer)

        return fprimer, rprimer if ret_str == False else str(fprimer), str(rprimer)





    @staticmethod
    def calcPhaseMatch(seq1, seq2, limit = 3):

        if not (isinstance(seq1, (str,Seq)) and isinstance(seq2, (str,Seq))):
            raise AssertionError('Seq 1 and Seq 2 must both be either a str or Seq object.')

        # if len(seq1) != len(seq2):
        #     raise AssertionError('Seq 1 and Seq 2 must be of the same length.')

        length = len(seq1) if len(seq1) == len(seq2) else max([len(seq1), len(seq2)])
        match_obj = SequenceMatcher(a=seq1, b=seq2)
        longest_match = match_obj.find_longest_match(0, length, 0, length)
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

        return True if phase_match_size <= limit else False






    @staticmethod
    def selfDimerizeTest(primer):
        if not isinstance(primer, Seq):
            primer = Seq(primer, generic_dna)
        length = len(primer)

        primer_Rev = primer[::-1]
        primer_Com = primer.complement()

        PLA_Seq.calcPhaseMatch(primer_Rev, primer_Com)




    @staticmethod
    def OrthogonalityTest(seq_1, seq_2, limit=4):
        if not (isinstance(seq_1, (str,Seq)) and isinstance(seq_2, (str,Seq))):
            raise AssertionError('Seq 1 and Seq 2 must both be either a str or Seq object.')

        #print(f"Checking orthogonality between:\n{seq_1} and \n{seq_2}")
        seq_1_Rev = seq_1[::-1]
        seq_2_Com = seq_2.complement()

        PLA_Seq.calcPhaseMatch(seq_1_Rev, seq_2_Com, limit)


    @staticmethod
    def orthogonalPoolComb(tuple_list, len_output_seq=0, round=1):
        print(250 * "_")
        all_combinations = combinations(tuple_list, 2)
        filtered_sequences = list()
        for pair1, pair2 in all_combinations:

            forwP_forwP = PLA_Seq.OrthogonalityTest(pair1.forw_primer, pair2.forw_primer, limit=5)
            forwP_revP = PLA_Seq.OrthogonalityTest(pair1.forw_primer, pair2.rev_primer, limit=5)
            revP_forwP = PLA_Seq.OrthogonalityTest(pair1.rev_primer, pair2.forw_primer, limit=5)
            revP_revP = PLA_Seq.OrthogonalityTest(pair1.rev_primer, pair2.rev_primer, limit=5)

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
            print(len_filtered_seq, round)
            return PLA_Seq.orthogonalPoolComb(filtered_sequences, len_output_seq=len_filtered_seq, round=round + 1)




def main():
    primer_pairs = PLA_Seq.genPrimerPairs()
    PLA_Seq.selfDimerizeTest(primer_pairs[0])
    # PLA_Seq.evalPrimerPairMT(*primer_pairs)
    #print(PLA_Seq.certMT_SelfDimerPrimerPairs())



# obj1 = PLA_Seq()
# obj2 = PLA_Seq(orth_seqs=['Hello'])
# obj3 = PLA_Seq()
#
# obj4 = PLA_Seq(orth_seqs=('Hello',))

#print(PLA_Seq.genPrimerPairs())

main()


