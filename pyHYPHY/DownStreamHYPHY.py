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

class WrongFileTypeForGapCheck(Exception):
    """
:exception
    :raise when given wrong file type
    """
    pass

class idSequenceUnKnow(Exception):
    """
    :raise when gene id sequence is fixed
    """
    pass

def get_gene_id_sequnce_from_lst(lst_file):
    """
    get gene_id from .lst file .
    .lst file was generated from grep procedure.
    So the sequence differs from gene_id generated from native Python os.listdir()
    :param lst_file:
    :return:
    """
    with open(lst_file, "r") as lst_fecther:
        line1 = lst_fecther.readline()
        lst_lines = lst_fecther.readlines()

    geneids = [line.split()[0] for line in lst_lines]
    return geneids


def get_gene_id(line, model_name):
    """
    extract gene id information from grep result of .result file.
    :param line:
    :param model_name:
    :return:
    """
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

def sum_available_sites(input_file_hyphy):

    if ".input" == os.path.splitext(input_file_hyphy)[-1]:
    # if input file is .input
        with open(input_file_hyphy, "r") as hyphy_input:
            genes = hyphy_input.readlines()[1::2]

        genes_no_tail = [gene.strip() for gene in genes]

    # if use .aln file
    elif ".aln" == os.path.splitext(input_file_hyphy)[-1]:
        with open(input_file_hyphy, "r") as reader:
            content = reader.readlines()

        genes_no_tail = [gene.split()[-1].strip() for gene in content]

    else:
        raise WrongFileTypeForGapCheck

    gene_ref = genes_no_tail[0]
    full_length = 0
    available_sites = 0

    for index, nt in enumerate(gene_ref):
        volume_matrix  = [gene[index] for gene in genes_no_tail]

        full_length += 1
        if not "-" in volume_matrix: # no gap in current site
            available_sites += 1

    return available_sites, full_length

def gap_check_traversal(input_folder, output_file, given_sequence_file="", jobids=[]):
    """
    check gap proportion on all files in "input_folder", and export as lst file to "output_file"
    with (optional) given_sequence_file

    :param input_folder:
    :param output_file:
    :param given_sequence:
    :return:
    """
    curdir_abs = os.path.abspath(os.curdir)

    if given_sequence_file:
        jobids = get_gene_id_sequnce_from_lst(given_sequence_file)
    elif jobids:
        pass
    else:
        raise idSequenceUnKnow
        #jobids = sorted([file_input for file_input in os.listdir(input_folder) if ".input" == os.path.splitext(file_input)[-1]])

    with open(output_file, "w") as writer:
        writer.write("geneid\tgaps\tfull\n")
        for jobid in jobids:
            input_name = os.path.join(input_folder, jobid + ".input")
            aln_name = os.path.join(input_folder, jobid + ".aln")
            if os.path.exists(input_name):
                num_gap, full_nt_length = gap_counting_input(input_name)
            elif os.path.exists(aln_name):
                num_gap, full_nt_length = gap_counting_input(aln_name)
            else:
                raise WrongFileTypeForGapCheck
            writer.write("%s\t%s\t%s\n" % (jobid, str(num_gap), str(full_nt_length)))

    os.chdir(curdir_abs)

def report_traversal_available_site(input_folder, output_file, given_sequence_file="", jobids=[]):
    """
    check available sites  proportion on all files in "input_folder", and export as lst file to "output_file"
    with (optional) given_sequence_file

    :param input_folder:
    :param output_file:
    :param given_sequence:
    :return:
    """
    curdir_abs = os.path.abspath(os.curdir)

    if given_sequence_file:
        jobids = get_gene_id_sequnce_from_lst(given_sequence_file)
    elif jobids:
        pass
    else:
        raise idSequenceUnKnow
        #jobids = sorted([file_input for file_input in os.listdir(input_folder) if ".input" == os.path.splitext(file_input)[-1]])

    with open(output_file, "w") as writer:
        writer.write("geneID\tavailableSites\tfullLength\n")
        for jobid in jobids:
            input_name = os.path.join(input_folder, jobid + ".input")
            aln_name = os.path.join(input_folder, jobid + ".aln")
            if os.path.exists(input_name):
                num_no_gap, full_nt_length = sum_available_sites(input_name)
            elif os.path.exists(aln_name):
                num_no_gap, full_nt_length = sum_available_sites(aln_name)
            else:
                raise WrongFileTypeForGapCheck
            writer.write("%s\t%s\t%s\n" % (jobid, str(num_no_gap), str(full_nt_length)))

    os.chdir(curdir_abs)

if __name__ == "__main__":
    pass