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


if __name__ == "__main__":
    main1()
    main2()
    main3()
    # num_gap, length_full = gap_counting_input("tmp.input")
    # print num_gap, length_full