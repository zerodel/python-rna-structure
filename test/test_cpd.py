# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#
#
import os
import os.path
import pyRNAsnp.pyRNAsnp as RNAsnp
import pyHYPHY.global_constants as gc
import pyHYPHY.ModelHYPHY as pyMH


import unittest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(False, False)

    def test_check_codons(self):
        ccl_file = "/home/zerodel/Workspace/codon_all.lst"
        mp = pyMH.ModelHYPHY(lst_2_tuple_file_name=ccl_file)
        feedback = mp.check_codons("AAA", "AAA")
        print str(feedback)
        self.assertTrue(tuple == type(feedback))



if __name__ == "__main__":
    unittest.main()