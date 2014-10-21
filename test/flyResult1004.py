# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#

import os
import os.path
import pyHYPHY.DownStreamHYPHY as dsh



dot_result_folder = "/home/zerodel/Workspace/Fly/FlyFull0930"

para_path = "/home/zerodel/Workspace/Fly/ExtractedParameters"

gtrmodel = "gtr"
nestmodel = "nest"

gtr_list_para = ["likelihood", "aic", "w"]
nest_list_para = ["likelihood", "aic", "w", "psi"]

def get_para_1004():
    # get parameters from fly data . model : gtr, nest, parameter: likelihood, aic, w, psi

    print "---- "
    dsh.extract_parameter(dot_result_folder, para_path, ["gtr"], gtr_list_para)
    dsh.extract_parameter(dot_result_folder, para_path, ["nest"], nest_list_para)

def prepare4r():
    dsh.transfer_for_r("gtr", "", gtr_list_para, "/home/zerodel/Workspace/Fly/ExtractedParameters")
    dsh.transfer_for_r("nest", "", nest_list_para, "/home/zerodel/Workspace/Fly/ExtractedParameters")


if __name__ == "__main__":
    get_para_1004()
    prepare4r()