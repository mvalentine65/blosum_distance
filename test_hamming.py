import unittest
from blosum_distance import cluster_distance_filter
from time import time

class TestClusterDistanceFilter(unittest.TestCase):
    def test_basic_functionality(self):
        lines = []
        dummy = 'dummy.fa'
        with open(dummy) as f:
            lines = f.readlines()
        print(f'Operating on {len(lines)} lines.')
        start = time()
        clusters = cluster_distance_filter(lines)
        end = time()
        print(f'Total time: {end-start}')
        

if __name__ == '__main__':
    unittest.main()
