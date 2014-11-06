# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
import os
import os.path
from math import log

__author__ = 'zerodel'

# 2014年10月29日， 现在mac 上git 下来的代码不是最新的，但是这些参数不变， 不影响操作。

models = ["gtr", "nest"]
np = {"gtr":7, "nest":8}
step_num = 20
step_width_step = 3
site_shift_step = 3

site_start_point = 0
window_width_min = 30

class SomeError(Exception):
    """ to indicate where an error occurs
"""
    pass


def export_parameter(result_folder, parameter_name, output_file, result_file_pattern="*.result"):

    current_dir = os.path.abspath(os.curdir)
    os.chdir(result_folder)
    ord_string = "grep -i %s %s > %s" % (parameter_name,result_file_pattern, output_file)
    print ord_string
    os.system(ord_string)
    os.chdir(current_dir)



def input_file_info(dot_input_file_name):
    """
    get info from .input file
"""
    length_seq = 0
    with open(dot_input_file_name) as reader:
        for line in reader:
            if not line.startswith(">"):
                len_this_seq = len(line)
                if 0 == length_seq:
                    length_seq = len_this_seq
                if not length_seq == len_this_seq:
                    # not same length in .input file
                    print "oops , not same length in .input file ...."
                    return None
    # if no error happens , it will return a interger showing the length of sequence in this .input file 
    return length_seq


def export_length(result_path,output_path):
    with open(os.path.join(output_path, "seqLen.lst"), "w") as seq_len_table:
        window_width_table_head = [str(step_width_step*single_offset + window_width_min)
                                   for single_offset in range(step_num)]
        seq_len_table.write("start_site/windowWidth\t"
                            + "\t".join(window_width_table_head) + "\n")
        for window_start_offset in range(step_num):
            start_site = site_shift_step*window_start_offset + site_start_point
            len_seq_input = []
            for window_width_offset in range(step_num):
                length_window = step_width_step*window_width_offset + window_width_min
                job_description = "w%ds%d" % (length_window, start_site)
                dot_input_file = os.path.join(result_path, "TIR%s.input" % job_description)
                length_seq_in_input = input_file_info(dot_input_file)
                if not isinstance(length_seq_in_input, int):
                    print "not same length , exiting..."
                    return
                len_seq_input.append(str(length_seq_in_input))
            seq_len_table.write("\t".join(len_seq_input) + "\n")


def clean_for_r_1(raw_file, output_path):
    """ data clean for HYPHY results , and prepared for R visualization
"""
    # do some data clean for R

    #output_file = "/Users/zerodel/WorkSpace/sWLikelihood.lst"
    with open(raw_file, "r") as result_reader:
        contents = result_reader.readlines()

    for modelname in models:
        writer_ML = open(os.path.join(output_path, "%sML.lst" % modelname), "w")
        writer_BIC= open(os.path.join(output_path, "%sBIC.lst" % modelname), "w")

        window_width_table_head = [str(step_width_step*single_offset + window_width_min)
                                   for single_offset in range(step_num)]

        table_head = "start_site/windowWidth\t" + "\t".join(window_width_table_head) + "\n"

        try:
            writer_ML.write(table_head)
            writer_BIC.write(table_head)

            for window_start_offset in range(step_num):
                ml_same_start_site = []
                bic_same_start_site = []
                start_site = site_shift_step*window_start_offset + site_start_point
                for window_width_offset in range(step_num):

                    length_window = step_width_step*window_width_offset + window_width_min

                    site_string = "s%d" % start_site
                    length_string = "w%d" % length_window
                    window_description = "w%ds%d" % (length_window, start_site)
                    job_description = "w%ds%d%s" % (length_window, start_site, modelname)


                    likelihood_needed = [line.strip().split("=")[-1] for line in contents if job_description in line]
                    if not len(likelihood_needed) == 1:
                        print site_string, "\t", length_string, "\t", modelname
                        print str(likelihood_needed)
                    # print likelihood_needed
                    likelihood_this_window = likelihood_needed[0].strip()
                    ml_same_start_site.append(likelihood_this_window)

                    # bic calculation
                    likelihood = float(likelihood_this_window)
                    np_model = np[modelname]
                    length_seq = input_file_info(os.path.join(result_path_nogap, "TIR%s.input" % window_description))
                    bic = np_model*log(length_seq) - 2.0*likelihood
                    bic_same_start_site.append(str(bic).strip())

                #print "\t".join(ml_same_start_site)
                writer_ML.write(str(start_site) + "\t" + "\t".join(ml_same_start_site) + "\n")
                writer_BIC.write(str(start_site) + "\t" + "\t".join(bic_same_start_site) + "\n")

        finally:
            writer_BIC.close()
            writer_ML.close()

def clean_for_r_other(raw_result_path, output_path, parameter_wanted, pattern_output="w%ds%d%s"):
    step_num = 21
    step_width_step = 15
    window_width_min = 30
    site_start_point = 0
    site_shift_step = 15

    window_axis = [0, 2]
    start_axis = range(step_num)

    lst_file_name_para = os.path.join(output_path, "%s.lst" % parameter_wanted)
    export_parameter(raw_result_path, parameter_wanted, lst_file_name_para)

    with open(lst_file_name_para, "r") as result_reader:
        contents = result_reader.readlines()

    for modelname in models:
        writer_ML = open(os.path.join(output_path, "%s%s.lst" % (modelname, parameter_wanted)), "w")

        window_width_table_head = [str(step_width_step*single_offset + window_width_min)
                                   for single_offset in window_axis]

        table_head = "start_site/windowWidth\t" + "\t".join(window_width_table_head) + "\n"

        try:
            writer_ML.write(table_head)


            for window_start_offset in start_axis:
                ml_same_start_site = []
                start_site = site_shift_step*window_start_offset + site_start_point
                for window_width_offset in window_axis:

                    length_window = step_width_step*window_width_offset + window_width_min

                    site_string = "s%d" % start_site
                    length_string = "w%d" % length_window
                    window_description = "w%ds%d" % (length_window, start_site)

                    job_description = pattern_output % (length_window, start_site, modelname)

                    # data grab
                    # para_needed = [line.strip().split("=")[-1] for line in contents if job_description in line]
                    para_needed = []
                    for line in contents:
                        if job_description in line and "=" in line:
                            para_needed.append(line.strip().split("=")[-1])

                        if job_description in line and "AIC" in line:
                            para_needed.append(line.strip().split()[-1])

                    if not len(para_needed) == 1:
                        print site_string, "\t", length_string, "\t", modelname
                        print str(para_needed)
                        if 0 == len(para_needed):
                            raise SomeError
                    # print likelihood_needed
                    para_this_window = para_needed[0].strip()
                    ml_same_start_site.append(para_this_window)
                writer_ML.write(str(start_site) + "\t" + "\t".join(ml_same_start_site) + "\n")
        except SomeError,e:
            print "model name: %s\t w:%d \t s:%d has no parameter named : %s" % (modelname, length_window, start_site, parameter_wanted)
            continue

        finally:
            writer_ML.close()

if __name__ == "__main__":
    result_path_nogap = "/Users/zerodel/WorkSpace/sliding300/nogap"
    parameter_wanted = "likelihood"
    output_folder = "/Users/zerodel/WorkSpace/sWAnalysis/300/nogap"
    output_file_likelihood = os.path.join(output_folder, "lnL.lst")
    #export_parameter(result_path_nogap, parameter_wanted, output_file_likelihood)
    clean_for_r_other(result_path_nogap, output_folder, "aic")
    #
    result_path_gap = "/Users/zerodel/WorkSpace/sliding300/gap"
    output_folder_gap = "/Users/zerodel/WorkSpace/sWAnalysis/300/gap"
    output_file_likelihood = os.path.join(output_folder_gap, "lnL.lst")
    # export_parameter(result_path, parameter_wanted, output_file_likelihood, "*N.result")
    clean_for_r_other(result_path_gap, output_folder_gap, "aic")
