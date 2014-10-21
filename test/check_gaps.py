# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#

import os
import os.path
import pyHYPHY.DownStreamHYPHY as dsh


def gap_counting_aln(alignment_file):
    with open(alignment_file, "r") as reader:
        content = reader.readlines()

    genes = [gene.split()[-1] for gene in content]
    full_nt = 0
    num_gap = 0
    for gene in genes:
        num_gap += gene.strip().count("-")
        full_nt += len(gene.strip())

    return num_gap, full_nt


def gap_counting_input(inputfile_hyphy):
    with open(inputfile_hyphy, "r") as reader:
        genes = reader.readlines()[1::2]
    full_nt = 0
    num_gap = 0
    for gene in genes:
        num_gap += gene.strip().count("-")
        full_nt += len(gene.strip())

    return num_gap, full_nt


def main1():
    ## get the gaps number and full_length .
    aln_path = "/home/zerodel/Workspace/Yeast/result/main_full_length"
    os.chdir(aln_path)
    aln_files = [file1 for file1 in os.listdir(aln_path) if ".aln" == os.path.splitext(file1)[-1]]
    jobids = dsh.get_gene_id_sequnce_from_lst(os.path.join("/home/zerodel/Workspace/Yeast/result/ExtractedParameter", "gtr.lst"))
    with open("/home/zerodel/Workspace/Yeast/result/ExtractedParameter/gapyeast.lst", "w") as writer:
        writer.write("gaps\tfull\n")
        for jobid in jobids:
            num_gap, full_nt_length = gap_counting_aln(jobid + ".aln")
            writer.write("%s\t%s\n" % (str(num_gap), str(full_nt_length)))


def main2():
    ## get the gaps number and full_length .
    aln_path = "/home/zerodel/Workspace/sp2"
    os.chdir(aln_path)
    aln_files = [file1 for file1 in os.listdir(aln_path) if ".aln" == os.path.splitext(file1)[-1]]
    jobids = dsh.get_gene_id_sequnce_from_lst(os.path.join("/home/zerodel/GitProjects/python-rna-structure/data/para", "nest2.lst"))
    with open("/home/zerodel/GitProjects/python-rna-structure/data/para/gap2.lst", "w") as writer:
        writer.write("gaps\tfull\n")
        for jobid in jobids:
            num_gap, full_nt_length = gap_counting_input(jobid + ".input")
            writer.write("%s\t%s\n" % (str(num_gap), str(full_nt_length)))


def main3():
    ## get the gaps number and full_length .
    aln_path = "/home/zerodel/Workspace/sp2"
    os.chdir(aln_path)
    aln_files = [file1 for file1 in os.listdir(aln_path) if ".aln" == os.path.splitext(file1)[-1]]
    jobids = dsh.get_gene_id_sequnce_from_lst(os.path.join("/home/zerodel/GitProjects/python-rna-structure/data/para", "nest2.lst"))
    with open("/home/zerodel/GitProjects/python-rna-structure/data/para/gap10.lst", "w") as writer:
        writer.write("gaps\tfull\n")
        for jobid in jobids:
            num_gap, full_nt_length = gap_counting_input(jobid + ".input")
            writer.write("%s\t%s\n" % (str(num_gap), str(full_nt_length)))


def gap_check_traversal(input_folder, output_file, given_sequence_file=""):
    """
    check gap proportion on all files in "input_folder", and export as lst file to "output_file"
    with (optional) given_sequence_file

    :param input_folder:
    :param output_file:
    :param given_sequence:
    :return:
    """
    curdir_abs = os.path.abspath(os.curdir)

    if not given_sequence_file:
        jobids = dsh.get_gene_id_sequnce_from_lst(os.path.join(given_sequence_file))
    else:
        raise dsh.idSequenceUnKnow
        #jobids = sorted([file_input for file_input in os.listdir(input_folder) if ".input" == os.path.splitext(file_input)[-1]])

    with open(output_file, "w") as writer:
        writer.write("gaps\tfull\n")
        for jobid in jobids:
            input_name = os.path.join(input_folder, jobid + ".input")
            aln_name = os.path.join(input_folder, jobid + ".aln")
            if os.path.exists(input_name):
                num_gap, full_nt_length = gap_counting_input(input_name)
            elif os.path.exists(aln_name):
                num_gap, full_nt_length = gap_counting_input(aln_name)
            else:
                raise dsh.WrongFileTypeForGapCheck
            writer.write("%s\t%s\n" % (str(num_gap), str(full_nt_length)))

    os.chdir(curdir_abs)

if __name__ == "__main__":
    # num_gap, length_full = gap_counting_input("tmp.input")
    # print num_gap, length_full
    pass