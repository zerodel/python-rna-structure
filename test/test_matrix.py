# !/usr/bin/env python
# -*- coding:utf-8 -*-
# Readme:
# changelog :nothing
#

__author__ = 'Zerodel_2'

import unittest
import pyHYPHY.global_constants as gc
import pyHYPHY.ModelHYPHY as pyMH
import pyHYPHY.MatrixHYPHY as pyMatrix
import os
import re

import logging
logging.basicConfig(filename="matrix_test.log", filemode="w", level=logging.DEBUG)

class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(False, False)

    def test_matrix(self):

        with open(name="gu_model.mdl", mode="r") as readerG:
            model_gu = [line for line in readerG if "{" in line or "}" in line]

        print len(model_gu)
        index_codon = gc.codon_list_hypothesis[:]

        mp = pyMH.ModelHYPHY()
        psi_in_model = 0
        for line_num, line in enumerate(model_gu):
            entry_in_line = line.split(",")
            for entry_index, single_entry in enumerate(entry_in_line):
                codon_check_response = mp.check_codons(codon_origin=index_codon[line_num], codon_target=index_codon[entry_index])

                if True == codon_check_response:
                    self.assertIn("*", single_entry)
                    logging.debug("\n codon1:%s \t codon2:%s \t entry:%s \n" % (index_codon[line_num], index_codon[entry_index], single_entry))
                elif type(codon_check_response) == int:
                    self.assertIn("0", single_entry)
                    logging.debug("\n codon1:%s \t codon2:%s \t entry:%s \n" % (index_codon[line_num], index_codon[entry_index], single_entry))

                elif type(codon_check_response) == tuple:
                    self.assertIn("a", single_entry)
                    logging.debug("\n codon1:%s \t codon2:%s \t entry:%s \n" % (index_codon[line_num], index_codon[entry_index], single_entry))

                    # 论证一下是否psi都是"*psi"形式存在
                    if "psi" in single_entry:
                        logging.debug("psi is in entry")
                        psi_in_model += 1
                        if "}" in single_entry:
                            pos = single_entry.index("}")
                            print single_entry, "\t", single_entry[0:pos], "\t", single_entry[pos:-1]
                    self.assertEqual("*psi" in single_entry, "psi" in single_entry)

            line_rebuilt = ",".join(entry_in_line)
            self.assertEqual(line, line_rebuilt)

        print " have %d psi in model" % psi_in_model

    def test_c2c(self):
        index = gc.codon_list_hypothesis
        mg = pyMH.ModelHYPHY(lst_2_tuple_file_name="../50CCL.lst")

        mg2 = pyMH.ModelHYPHY()
        tuple_list = mg2.read_lst("../50CCl.lst")

        all_CCLs = 0
        for i1, c1 in enumerate(index):
            for i2, c2 in enumerate(index):
                response = mg.check_codons(c1, c2)
                if tuple == type(response):
                    if response[2]:
                        all_CCLs += 1
                        self.assertIn((c1, c2), tuple_list)
                        logging.debug("c1:%s \t c2:%s" % (c1, c2))

        self.assertEqual(all_CCLs, len(tuple_list))
        print len(tuple_list), all_CCLs

    def test_nest_export(self):
        pyMatrix.degenerate(gtr_nested_file="gu_model.mdl", gtr_mdl_file="rebuilt_gtr.mdl")
        model_your = pyMatrix.nest_export(gtr_mdl="rebuilt_gtr.mdl", gtr_nested="rna_structure.mdl", ccl_file="../50CCl.lst")
        print len(re.findall("\*psi", model_your)), "psi in model_your"
        self.assertTrue(os.path.exists("rna_structure.mdl"))

        index_codon = gc.codon_list_hypothesis
        mg2 = pyMH.ModelHYPHY()
        tuple_list = mg2.read_lst("../50CCl.lst")

        zero_psi = True
        with open(name="rna_structure.mdl", mode="r") as readerG:
            model_nest = [line for line in readerG ]

        print len(model_nest), " lines in model"
        for line_num, line in enumerate(model_nest):
            entry_in_line = line.split(",")
            for entry_index, single_entry in enumerate(entry_in_line):
                if "psi" in single_entry:
                    zero_psi = False
                    self.assertIn((index_codon[line_num], index_codon[entry_index]), tuple_list)

        self.assertFalse(zero_psi)


if __name__ == '__main__':
    unittest.main()