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

def get_gene_id_sequnce_from_lst(lst_file):
    with open(lst_file, "r") as lst_fecther:
        line1 =lst_fecther.readline()
        lst_lines = lst_fecther.readlines()

    geneids = [line.split()[0] for line in lst_lines]
    return geneids

def get_gene_id(line, model_name):
    return line.split(":")[0].split(".")[0][: - len(model_name)]


def file_extract(file_name, model_name):
    """
    read from para/*.txt and return gene_id and parameter value from string like
    "b4200nest.result:w=0.173206"
    :param file_name:
    :return: two list : gene_ids, para_values
    """
    with open(file_name, "r") as reader:
        contents = reader.readlines()

    gene_ids = [get_gene_id(line, model_name) for line in contents]
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
    para_name = ["geneid"]
    for para in para_list:
        file_name_tmp = os.path.join(folder_path, model_name + mid + para + ".txt")
        gene_index, para_content = file_extract(file_name_tmp, model_name)

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


def extract_parameter(result_path, para_path, list_model, list_parameter, mid=""):
    #use grep to get parameters from .result files
    dir_now = os.path.abspath(os.curdir)

    os.chdir(result_path)

    for para in list_parameter:
            for model in list_model:
                system_order = "grep -i " + para + " *" + model + ".result >" + os.path.join(para_path, "".join([model, mid, para, ".txt"]))
                os.system(system_order)

    os.chdir(dir_now)


def test_gene_species_match(aln_file_path, list_gene_given):
    aln_files = [single_file for single_file in os.listdir(aln_file_path) if ".aln" == single_file[-4:]]
    num_not_match = 0
    genes_matched = []
    print "there are ", str(len(aln_files)), "genes here"
    for file_aln_1 in aln_files:
        full_path = os.path.join(aln_file_path, file_aln_1)
        with open(full_path, "r") as reader:
            contents = reader.readlines()
            species = [line.split()[0].strip() for line in contents]
        match_newick_species = True
        for single_species in species:
            if not single_species in list_gene_given:
                print "some gene not match in ", file_aln_1
                match_newick_species = False

        if match_newick_species:
            print file_aln_1, "--- matched !"
            genes_matched.append(file_aln_1.split(".")[0])
        else:
            num_not_match += 1

    print "In summary , ", str(num_not_match), "genes unmatch"


if __name__ == "__main__":
    pass