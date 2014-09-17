# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#

import pyHYPHY.DownStreamHYPHY as dsh
import os
import os.path

def main2():
    path_result = "/home/zerodel/GitProjects/python-rna-structure/data/para"
    path_output = "/home/zerodel/Workspace/former10"

    dir_now = os.path.abspath(os.curdir)

    os.chdir(path_result)

    paras = ["likelihood", "aic", "w"]
    mid = "10"
    models = ["gtr", "gu", "nest"]

    for para in paras:
            for model in models:
                system_order = "grep -i " + para + " *" + model + ".result >" + "../" + model + mid + para + ".txt"
                os.system(system_order)

    os.chdir(dir_now)


def main3():
    path_result = "/home/zerodel/Workspace/merge0914/2sp"

    dir_now = os.path.abspath(os.curdir)

    os.chdir(path_result)
    para = "psi"
    model = "nest"
    mid = "2"
    system_order = "grep -i " + para + " *" + model + ".result >" + "../" + model + mid + para + ".txt"
    os.system(system_order)
    path_working = "/home/zerodel/Workspace/merge0914"

    gene_all = dsh.transfer_for_r("nest", mid, ["psi"], path_working)


def main1():
    paras = ["likelihood", "aic", "w"]
    path_working = "/home/zerodel/Workspace/former10"
    mid = "10"

    gene_all = dsh.transfer_for_r("gtr", mid, paras, path_working)
    gene_all_tmp = dsh.transfer_for_r("gu", mid, paras, path_working)

    if not gene_all == gene_all_tmp:
        print "not same length"

    gene_all_tmp = dsh.transfer_for_r("nest", mid, paras, path_working)

    if not gene_all == gene_all_tmp:
        print "not same length"

    #
    # mid = "2"
    # gene_all = dsh.transfer_for_r("gtr", mid, paras, path_working)
    #
    # gene_all_tmp = dsh.transfer_for_r("gu", mid, paras, path_working)
    # if not gene_all == gene_all_tmp:
    #     print "not same length"
    # gene_all_tmp = dsh.transfer_for_r("nest", mid, paras, path_working)
    # if not gene_all == gene_all_tmp:
    #     print "not same length"


if __name__ == "__main__":
    main3()