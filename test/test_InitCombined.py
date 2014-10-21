__author__ = 'Zerodel'

import unittest

from InitCombined import *
import pyHYPHY.DataHYPHY as DH
import os
import os.path
class MyTestCase(unittest.TestCase):
    def test_paste_matrix(self):
        m1 = ["123", "abc"]
        m2 = ["x", "0"]
        mm = ["123x", "abc0"]
        self.assertEqual(mm, paste_matrix(m1,m2))


    def test_extract_TIR(self):
        filename = "tmp.aln"
        TIR_len = 20

        matrix_TIR = extract_TIR_single_file(filename, TIR_len)
        print matrix_TIR

        mm = ["12345678901234567", "12345678901234456"]
        mr = DH.remove_gaps_matrix(matrix_TIR)
        print mr
        self.assertLess(len(mr[0]), len(mm[0]))


    def test_gene_conjunction(self):
        path_aln = "."
        inputfile = "combinedTIR.input"
        lengthTIR = 20
        gene_conjunction(path_aln,inputfile,lengthTIR)

        self.assertTrue(os.path.exists(inputfile))




if __name__ == '__main__':
    unittest.main()
