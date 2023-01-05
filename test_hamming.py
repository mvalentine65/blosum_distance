import os
import unittest

import numpy as np
import phymmr_tools as phy
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pprint import pprint
from time import time
from itertools import combinations
from random import randint

def make_indices(sequence: str, gap_character="-") -> tuple:
    """
    Finds the index of the first and last non-gap bp in a sequence.
    Returns the start value and the end values + 1 as a tuple.
    """
    start = None
    end = None
    for i, character in enumerate(sequence):
        if character != gap_character:
            start = i
            break
    for i in range(len(sequence) - 1, -1, -1):
        if sequence[i] != gap_character:
            end = i + 1
            break
    if start is None or end is None:
        raise ValueError()
    return start, end


def find_distance_from_matrix(thing, start, end):
    cross, maxtuple = thing
    try:
        if start != 0:
            return 1 - ((cross[end - 1] - cross[start - 1]) / max(maxtuple[0][end - 1] - maxtuple[0][start - 1],
                                                maxtuple[1][end - 1] - maxtuple[1][start - 1]))
        else:
            return 1 - (cross[end - 1] / max(maxtuple[0][end - 1], maxtuple[1][end - 1]))
    except ZeroDivisionError:
        return np.nan

class TestReferenceMatrix(unittest.TestCase):

    def setUp(self):
        # self.gene = "600.fa"
        self.gene ="JSM_2594_S5.fa/trimmed/aa"
        # self.files = "JSM_2594_S5.fa/trimmed/aa/EOG7B0GZN.aa.fa"
        self.files = [os.path.join(self.gene, x) for x in os.listdir(self.gene)]

    def test_reference_matrix_core_function(self):
        comb = 0
        py = 0
        start = 5
        end = 50
        for inf in self.files:
            with open(inf) as f:
                seq_list = [x for x in SimpleFastaParser(f)]
            headers = [x[0] for x in seq_list if x[0][-1] == "."]
            ref_seqs = [x[1] for x in seq_list if x[0][-1] == "."]
            matrix_start = time()
            ref_distance_matrix = phy.make_ref_distance_matrix_supervector(ref_seqs)
            matrix_time = time() - matrix_start
            length = len(ref_seqs[0])
            man_start = time()
            for i in range(10):
                manual = []
                abridged = [x[1][start:end] for x in seq_list]
                for i in range(len(headers)-1):
                    for j in range(i+1,len(headers)):
                        manual.append(phy.blosum62_distance(abridged[i], abridged[j]))
            man_end = time()
            py += man_end-man_start
            inpy_start = time()
            for i in range(10):
                start = 5
                end = 50
                output = [find_distance_from_matrix(seq, start, end) for seq in ref_distance_matrix]
            inpy_end = time()
            comb += inpy_end-inpy_start
            manual2 = [x for x in manual if not np.isnan(x)]
            output2 = [x for x in output if not np.isnan(x)]
            if manual2 and output2:
                q1_man = np.nanpercentile(manual2, 25, method="midpoint")
                q3_man = np.nanpercentile(manual2, 75, method="midpoint")

                q1_comb = np.nanpercentile(output2, 25, method="midpoint")
                q3_comb = np.nanpercentile(output2, 75, method="midpoint")
                pass
                if q1_man != q1_comb:
                    with open(f'{start}-{end}.csv','a+') as out:
                        out.write(f'{inf},{q1_man},{q1_comb}\n{str(manual2)}\n{str(output2)}\n\n')
                # self.assertAlmostEqual(q1_man, q1_comb,2)
            # self.assertEqual(manual, subtracted)
            # self.assertEqual((q1_man, q3_man), (q1_sub, q3_sub))

        print(f'com_time: {comb}')
        print(f'man_time: {py}')
        print("\n")


class TestIndices(unittest.TestCase):

    def setUp(self):
        self.file = "JSM_2594_S5.fa/trimmed/aa/EOG7B0H10.aa.fa"

    def test_indices(self):
        with open(self.file) as f:
            seq_list = [x for x in SimpleFastaParser(f)]
        seqs = [x[1] for x in seq_list]
        reps = 10000
        py_start = time()
        for i in range(reps):
            for seq in seqs:
                make_indices(seq)
        py_end = time()
        rs_start = time()
        for i in range(reps):
            for seq in seqs:
                phy.make_indices(seq)
        rs_end = time()
        print(f'py: {py_end-py_start}')
        print(f'rs: {rs_end-rs_start}')

    def test_index_correctness(self):
        with open(self.file) as f:
            seq_list = [x for x in SimpleFastaParser(f)]
        seqs = [x[1] for x in seq_list]
        for seq in seq_list:
            self.assertEqual(make_indices(seq), phy.make_indices(seq))



if __name__ == '__main__':
    unittest.main()
