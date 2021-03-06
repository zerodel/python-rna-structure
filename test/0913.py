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


class NotConsistentLength(Exception):
    """
    exception thrown  when the entry length are different
    """
    pass


def file_extract(file_name):
    """
    read from para/*.txt and return gene_id and parameter value from string like
    "b4200nest.result:w=0.173206"
    :param file_name:
    :return: two list : gene_ids, para_values
    """
    with open(file_name, "r") as reader:
        contents = reader.readlines()

    gene_ids = [line[:5] for line in contents]
    if "=" in contents[-1]:
        para_values = [line.split("=")[-1] for line in contents]
    else:
        # aic value is different from other parameter, my fault
        para_values = [line.split()[-1] for line in contents]
    return gene_ids, para_values


def data_grab(model_name, mid, para_list, folder_path):
    """
    grab data from .txt files
    :param model_name: like "gtr", "gu", "nest"
    :param mid: 10 or 2
    :param para_list: "w" "aic" "likelihood"
    :param folder_path: path to .txt files
    :return:
    """
    gene_index_matrix = []
    para_matrix = []
    para_name = ["gene_id"]
    for para in para_list:
        file_name_tmp = os.path.join(folder_path, model_name + mid + para + ".txt")
        gene_index, para_content = file_extract(file_name_tmp)

        para_name.append(para)
        gene_index_matrix.append(gene_index)
        para_matrix.append(para_content)

    gene_index = gene_index_matrix[0]
    for gene_index_given in gene_index_matrix:
        if not gene_index == gene_index_given:
            raise NotConsistentLength

    return gene_index, para_name, para_matrix


def transfer_for_r(modelname, mid, paras, path_working):
    """
    write .lst file , convenient for R program to import . arguments are the same as "data_grab"
    :param modelname:
    :param mid:
    :param paras:
    :param path_working:
    :return:
    """
    try:
        gene_index, para_name, para_matrix = data_grab(modelname, mid, paras, path_working)

    except NotConsistentLength:
        print "sorry , something wrong , num of entry not same length "
        return

    filename_writen_to = os.path.join(path_working, modelname + mid + ".lst")
    with open(filename_writen_to, 'w') as writer:
        writer.write("\t".join(para_name) + "\n")
        for index_gene, gene_single in enumerate(gene_index):
            # assemble single data line for data.frame
            tmp_data = [paras[index_gene].strip() for paras in para_matrix]
            data_line = gene_single + "\t" + "\t".join(tmp_data) + "\n"
            writer.write(data_line)

    return gene_index

def main1():
    paras = ["likelihood", "aic", "w"]
    mid = "10"
    modelgtr = "gtr"
    path_working = "/home/zerodel/GitProjects/python-rna-structure/data/para"
    gene_all = transfer_for_r("gtr", mid, paras, path_working)

    gene_all_tmp = transfer_for_r("gu", mid, paras, path_working)

    if not gene_all == gene_all_tmp:
        print "not same length"

    gene_all_tmp = transfer_for_r("nest", mid, paras, path_working)

    if not gene_all == gene_all_tmp:
        print "not same length"


    mid = "2"
    gene_all = transfer_for_r("gtr", mid, paras, path_working)

    gene_all_tmp = transfer_for_r("gu", mid, paras, path_working)
    if not gene_all == gene_all_tmp:
        print "not same length"
    gene_all_tmp = transfer_for_r("nest", mid, paras, path_working)
    if not gene_all == gene_all_tmp:
        print "not same length"


def main2():
    cpd_main = RNAsnp.get_whole_cpd("/home/zerodel/Workspace/cpd_store_site")
    codon_vector = gc.codon_list_hypothesis
    with open("/home/zerodel/GitProjects/python-rna-structure/data/rnasnp.matrix", "w") as writer:
        for codon1 in codon_vector:
            # build a line for lst file
            line_single = []
            for codon2 in codon_vector:
                line_single.append(str(len(cpd_main[(codon1, codon2)])))

            writer.write("\t".join(line_single) + "\n")
    with open("/home/zerodel/GitProjects/python-rna-structure/data/codon_vector.txt", "w") as writer:
        writer.write("\n".join(gc.codon_list_hypothesis))


def main3():
    cpd_main = RNAsnp.get_whole_cpd("/home/zerodel/Workspace/cpd_store_site")
    codon_vector = gc.codon_list_64
    with open("/home/zerodel/GitProjects/python-rna-structure/data/fullRnasnp.matrix", "w") as writer:
        available_cc = cpd_main.keys()
        for codon1 in codon_vector:
            # build a line for lst file
            line_single = []
            for codon2 in codon_vector:
                if (codon1, codon2) in available_cc:
                    line_single.append(str(len(cpd_main[(codon1, codon2)])))
                else:
                    line_single.append(str(0))

            writer.write("\t".join(line_single) + "\n")
    with open("/home/zerodel/GitProjects/python-rna-structure/data/codon_vector64.txt", "w") as writer:
        writer.write("\n".join(gc.codon_list_64))


if __name__ == "__main__":
    main3()