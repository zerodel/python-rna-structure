__author__ = 'Zerodel'

import unittest
import pyHYPHY.DataHYPHY as dh

class MyTestCase(unittest.TestCase):
    def test_something(self):
        raw_sequence = "aca----gt"
        gap_removed = "aca"

        self.assertEqual(gap_removed, dh.remove_gaps(raw_sequence))

    def test_remove_gap_matrix(self):
        raw_matrix = ["aaa--a", "bbbcac"]
        fileterd_matrix = ["aaa", "bbb"]

        self.assertEqual(fileterd_matrix, dh.remove_gaps_matrix(raw_matrix))

    def test_aln_to_no_gap_input(self):
        try:
            dh.check_sequence_matrix(["aaca", "aaa"])
        except dh.sequenceNotSameLength:
            pass
        else:
            self.fail("here should raise a exception!")

        dh.aln2inputNogap("./tmp.aln","./tmp2.input")




if __name__ == '__main__':
    unittest.main()
