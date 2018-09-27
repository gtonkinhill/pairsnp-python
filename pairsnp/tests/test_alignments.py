import unittest
from pairsnp import *
import numpy as np

class TestPairsnp(unittest.TestCase):
    def test_ambig(self):
        sparse_matrix, consensus, seq_names = calculate_snp_matrix("./pairsnp/tests/ambig.aln")
        d = calculate_distance_matrix(sparse_matrix, consensus, "dist", False)
        self.assertTrue(np.array_equal(d, np.array([[0,   0,  2,  1,  1],
                                        [0, 0,  2,  2,  2],
                                        [2, 2,  0,  3,  3],
                                        [1, 2,  3,  0,  0],
                                        [1, 2,  3,  0,  0]])))

    def test_ambig_with_n(self):
        sparse_matrix, consensus, seq_names = calculate_snp_matrix("./pairsnp/tests/ambig.aln")
        d = calculate_distance_matrix(sparse_matrix, consensus, "dist", True)
        self.assertTrue(np.array_equal(d, np.array([[0, 2,  4,  3,  3],
                                                    [2, 0,  4,  4,  4],
                                                    [4, 4,  0,  5,  5],
                                                    [3, 4,  5,  0,  0],
                                                    [3, 4,  5,  0,  0]])))

    def test_empty(self):
        self.failUnlessRaises(ValueError, lambda:calculate_snp_matrix("./pairsnp/tests/empty.aln"))

    def test_singleton(self):
        sparse_matrix, consensus, seq_names = calculate_snp_matrix("./pairsnp/tests/singleton.aln")
        d = calculate_distance_matrix(sparse_matrix, consensus, "dist", True)
        self.assertTrue(np.array_equal(d, np.array([[0]])))

    def test_bad(self):
        self.failUnlessRaises(ValueError, lambda:calculate_snp_matrix("./pairsnp/tests/bad.aln"))

    def test_lowercase(self):
        sparse_matrix, consensus, seq_names = calculate_snp_matrix("./pairsnp/tests/lowercase.aln")
        d = calculate_distance_matrix(sparse_matrix, consensus, "dist", False)
        self.assertTrue(np.array_equal(d, np.array([[0, 1,  2,  3],
                                                    [1, 0,  3,  4],
                                                    [2, 3,  0,  4],
                                                    [3, 4,  4,  0]])))

