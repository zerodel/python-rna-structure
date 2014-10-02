# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#

import os
import os.path
import pyHYPHY.DownStreamHYPHY as dsh

data_path = "/home/zerodel/Workspace/Yeast/result/main_full_length"
para_path = "/home/zerodel/Workspace/Yeast/result/ExtractedParameter"

def check_ccl():
    ccl1 = "/home/zerodel/Workspace/codon_all.lst"
    ccl2 = "/home/zerodel/Workspace/Yeast/codon_all.lst"
    with open(ccl1, "r") as cclreader:
        ccl1s = cclreader.read().split()

    with open(ccl2, "r") as cclreader:
        ccl2s = cclreader.read().split()

    print len(ccl1s), len(ccl2s)



if __name__ == "__main__":
    # dsh.extract_parameter(data_path, para_path, ["nest"], ["likelihood", "aic", "w", "psi"])
    # dsh.extract_parameter(data_path, para_path, ["gtr"], ["likelihood", "aic", "w"])
    # 
    # dsh.transfer_for_r("gtr", "", ["likelihood", "aic", "w"], para_path)
    # dsh.transfer_for_r("nest", "", ["likelihood", "aic", "w", "psi"], para_path)
    check_ccl()