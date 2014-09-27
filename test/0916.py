# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#

import pyRNAsnp.pyRNAsnp as RNAsnp
import pyHYPHY.DataHYPHY as DataHyPHY
import os
import os.path
import global_constants as gc
import pyHYPHY.MatrixHYPHY as matrixh

import pickle

__author__ = 'zerodel'


def test_cpd():
    cpd_path_name = "/home/zerodel/Workspace/cpd_store_site"

    geneid = "b0002"
    full_name = os.path.join(cpd_path_name, geneid + ".cpd")
    with open(full_name, "r") as cpdf:
        cpd_temp = pickle.load(cpdf)

    for cpdkey in cpd_temp.keys():
        print len(cpd_temp[cpdkey])


if __name__ == "__main__":
    test_cpd()