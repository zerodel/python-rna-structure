__author__ = 'Zerodel'

import unittest
import os
import os.path
import pyHYPHY.DownStreamHYPHY as dsh



class MyTestCase(unittest.TestCase):
    def test_gap(self):
        file_hyphy_input = "./tmp.aln"
        availble_sites,full_length = dsh.sum_available_sites(file_hyphy_input)
        self.assertEqual(full_length, 4)
        self.assertEqual(availble_sites,3)


    def test_gap_check_traversal(self):
        folder_input = os.path.abspath(os.curdir)
        output_gap_lst = "gap_report.txt"
        try:
            dsh.gap_check_traversal(folder_input, output_gap_lst)
        except dsh.idSequenceUnKnow:
            print "ok , exception raisen"
        else:
            self.fail("nope , no exception when there should have some ")

    def test_gap_traversal2(self):
        folder_input = os.path.abspath(os.curdir)
        output_gap_lst = "gap_report.txt"
        jobids = ["tmp"]
        dsh.gap_check_traversal(input_folder=folder_input, output_file=output_gap_lst,jobids=jobids)
        self.assertTrue(os.path.exists(output_gap_lst))

    def test_available_traversal(self):
        folder_input = os.path.abspath(os.curdir)
        output_gap_lst = "available_sites.txt"
        job_ids = ["tmp"]
        try:
            dsh.report_traversal_available_site(input_folder=folder_input,output_file=output_gap_lst)
        except dsh.idSequenceUnKnow:
            pass
        else:
            self.fail("id unknow")

        dsh.report_traversal_available_site(input_folder=folder_input,output_file=output_gap_lst,
        jobids=job_ids)

        self.assertTrue(os.path.exists(output_gap_lst))


if __name__ == '__main__':
    unittest.main()
