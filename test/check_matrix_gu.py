# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#
__author__ = 'Zerodel'



def AddSysPath(new_path):
    """ AddSysPath(new_path): adds a "directory" to Python's sys.path
    Does not add the directory if it does not exist or if it's already on
    sys.path. Returns 1 if OK, -1 if new_path does not exist, 0 if it was
    already on sys.path.
    """
    import sys, os
    # Avoid adding nonexistent paths
    if not os.path.exists(new_path): return -1
    # Standardize the path.  Windows is case-insensitive, so lowercase
    # for definiteness if we are on Windows.
    new_path = os.path.abspath(new_path)
    if sys.platform == 'win32':
        new_path = new_path.lower( )
        # Check against all currently available paths
        for x in sys.path:
            x = os.path.abspath(x)
            if sys.platform == 'win32':
                x = x.lower( )
                if new_path in (x, x + os.sep):
                    return 0
                sys.path.append(new_path)
                # if you want the new_path to take precedence over existing
                # directories already in sys.path, instead of appending, use:
                # sys.path.insert(0, new_path)
                return 1

AddSysPath("..")


try:
    import pyHYPHY.global_constants as hyphygc
    import pyHYPHY.ModelHYPHY as MH
except ImportError:
    print "still, can not import  "

# set the index of matrix , as alphabet

import unittest
import copy

index_hypothesis = copy.copy(hyphygc.codon_list_hypothesis)
aa_codon_dict = copy.copy(hyphygc.codon_aa_dict)

line1_matrix_gu = "*,aAC*w,aAG*psi,aAT*w,aAC*w,0,0,0,aAG*w,0,0,0,aAT*w,0,0,0,aAC*w,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,aAG*w,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"

class MyTestCase(unittest.TestCase):
    def test_matrix_index_line1(self):
        model_checker = MH.ModelHYPHY()
        print "start check"

        read_gu_model = open("test.mdl")
        lines = read_gu_model.readlines()

        for row_index, line1_matrix_gu in enumerate(lines):
            line1_matrix_gu_cleaned = line1_matrix_gu.translate(None,"{};")
            line1_entrys = [line_entry.strip() for line_entry in line1_matrix_gu_cleaned.split(",")]

            for index_column, entry_matrix in enumerate(line1_entrys):
                differ_codon = model_checker.check_codons(index_hypothesis[row_index], index_hypothesis[index_column])

                # print index_hypothesis[0], index_hypothesis[index_column], "differs ---> ", line1_entrys[index_column]
                # print str(differ_codon)
                if isinstance(differ_codon, bool):
                    self.assertEqual("*", line1_entrys[index_column])
                    continue
                if isinstance(differ_codon, int):  # differ more than one position
                    self.assertEqual("0", line1_entrys[index_column])

                if isinstance(differ_codon, tuple):
                    is_trans, is_sys, is_larger, nt_diff = differ_codon
                    self.assertIn(nt_diff, line1_entrys[index_column])
                    if not aa_codon_dict[index_hypothesis[row_index]] == aa_codon_dict[index_hypothesis[index_column]]:
                        self.assertIn("w", line1_entrys[index_column])

        read_gu_model.close()



if __name__ == "__main__":
    unittest.main()