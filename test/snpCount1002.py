# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#
import os
import os.path
import pyHYPHY.DownStreamHYPHY as dsh


header_mode3 = "SNP     w       Slen    GC      interval        d       pvalue1 ewin    interval        d_max   pvalue2"


# dsh.get_gene_id_sequnce_from_lst()

def count_snp_per_file(dot_rnasnp, is_greater=False, threshold=0.1, criterion="pvalue2"):
    # read a .rnasnp file and return snp count for this gene .
    # see Gu2014 for more information (We define structurally 88 sensitive sites in mRNA as those with putative local 89 structure-disruptive mutations)

    # 1st, read .rnasnp file
    with open(dot_rnasnp, "r") as rna_fetch:
        header = rna_fetch.readline().split()
        lines = rna_fetch.readlines()

    position_criterion = header.index(criterion)
    data_criterion = [float(line.split()[position_criterion]) for line in lines]

    if is_greater:
        result = [data_entry for data_entry in data_criterion if data_entry > threshold]
    else:
        result = [data_entry for data_entry in data_criterion if data_entry < threshold]

    return len(result)


def snp_count_species(lst_file, folder_rnasnp_file):
    # read gene id sequence
    path_dir = os.path.dirname(os.path.abspath(lst_file))
    lst_name_bare = os.path.splitext(lst_file)[0]

    jobids = dsh.get_gene_id_sequnce_from_lst(lst_file)
    file_name_vector = [os.path.join(folder_rnasnp_file, ".".join([jobid, "rnasnp"])) for jobid in jobids]

    is_greater_than = False
    data_criterion = "pvalue2"
    threshold = 0.1

    gene_snp_count = []
    for single_file in file_name_vector:
        if os.path.exists(single_file):
            gene_snp_count.append(str(count_snp_per_file(single_file, is_greater_than, threshold, data_criterion)/3.0))
        else:
            gene_snp_count.append("NA")
    # writer a new lst file for snp count in some species
    with open(os.path.join(path_dir, lst_name_bare + "SnpCount.lst"), "w") as exporter:
        exporter.write("geneid\tsnpCountPerGene\n")
        for index_g, jobid in enumerate(jobids):
            line_to_export = "%s\t%s\n" % (jobid, gene_snp_count[index_g])
            exporter.write(line_to_export)


if __name__ == "__main__":
    snp_count_species("/home/zerodel/Workspace/Yeast/result/ExtractedParameter/gtr.lst",
                      "/home/zerodel/Workspace/Yeast/yeast_rnasnp")

    snp_count_species("/home/zerodel/GitProjects/python-rna-structure/data/para/gtr10.lst",
                      "/home/zerodel/Workspace/ecoli_snp")