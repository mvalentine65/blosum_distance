import unittest
from blosum_distance import seqs_within_distance
from blosum_distance import str_hamming
from time import time

class TestDistanceCorrectness(unittest.TestCase):

    def test_seqs_within_no_diff_passes(self):
        one = "ATCG"
        two = "ATCG"
        self.assertTrue(seqs_within_distance(one, two, 1))

    def test_seqs_within_1_off_passes(self):
        one = "ATCG"
        two = "ATGG"
        self.assertTrue(seqs_within_distance(one, two, 1))

    def test_seqs_within_2_off_fails(self):
        one = "ATCG"
        two = "ATGC"
        self.assertFalse(seqs_within_distance(one, two, 1))

    def test_seqs_within_length_diff_fails(self):
        one = "MVEMJSUNP"
        two = "MVEMJSUN" # heresy
        self.assertFalse(seqs_within_distance(one, two, 1))

    def test_str_hamming_no_diff_passes(self):
        one = "ATCG"
        two = "ATCG"
        self.assertTrue(str_hamming(one, two, 1))

    def test_str_hamming_1_off_passes(self):
        one = "ATCG"
        two = "AACG"
        self.assertTrue(str_hamming(one, two, 1))

    def test_str_hamming_2_off_fails(self):
        one = "ATCG"
        two = "ACAG"
        self.assertFalse(str_hamming(one, two, 1))

    def test_str_hamming_length_diff_fails(self):
        one = "ATCG"
        two = "ACAGG"
        self.assertFalse(str_hamming(one, two, 1))

        
class TestDistanceBenchmarks(unittest.TestCase):

    def test_speed_identical_seqs(self):
        test_seq = "A"*50000
        reps = 10000
        seqs_start = time()
        for i in range(reps):
            _ = seqs_within_distance(test_seq, test_seq, 1)
        seqs_end = time()
        print(f"seqs_within_distance time: {seqs_end - seqs_start}")
        str_start = time()
        for i in range(reps):
            _ = str_hamming(test_seq, test_seq, 1)
        str_end = time()
        print(f"str_hamming time: {str_end - str_start}")


if __name__ == "__main__":
    unittest.main()
