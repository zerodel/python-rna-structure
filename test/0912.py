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


def check_aln_files():
    """
    check which .aln file contains gene of all 10 species
    :return:
    """
    folder_alns = "/media/zerodel/Home/Work/custom/ecoli_aln"
    all_alns = [single_file for single_file in os.listdir(folder_alns) if single_file[-4:] == ".aln"]
    aln_gene_num = []
    for single_file in all_alns:
        full_aln_path = os.path.join(folder_alns, single_file)
        genes, gene_lengths = DataHyPHY.aln_info(full_aln_path)
        aln_gene_num.append(len(genes))

    max_num = max(aln_gene_num)
    aln_full_gene = []
    for indexI, aln in enumerate(all_alns):
        if aln_gene_num[indexI] == max(aln_gene_num):
            aln_full_gene.append(all_alns[indexI])

    print "fulllength has", len(aln_full_gene), "with " , str(max_num), "genes"


def test_group_divide():
    path = "/home/zerodel/Workspace/codon_lst"
    divide_file = "/home/zerodel/Workspace/codon_all.lst"
    RNAsnp.group_divide(path, divide_file)




def built_model():
    parent_path = "/home/zerodel/GitProject"
    work_path = "/home/zerodel/Workspace"
    source_path = "/home/zerodel/Workspace/ecoli_snp"
    output_path = "/home/zerodel/Workspace/cpd_store_site"
    single_codon_file_path = os.path.join(work_path, "codon_lst")
    cc_significant_file_path = os.path.join(work_path, "codon_all.lst")
    rnasnp_files, snppathname = RNAsnp.snp_dir_list(source_path)



    #step 1 , set up  .cpd files
    # RNAsnp.snp_dir_traversal(rnasnp_files, snppathname, output_path, True)
    # # step 2 , print each codon lst file
    # for codon in gc.codon_list_hypothesis:
    #     codon_with_direction = codon + "_"
    #     print "--->", codon_with_direction
    #     RNAsnp.cpd_dir_traversal(codon_with_direction, output_path, single_codon_file_path, 60)



    # kmeans meet some problem here 20140912
    # # R kmeans part ******
    # print "kmeans of cluster"
    # kmeans_script_path = os.path.join(parent_path,"kmeans.R")
    # cc_significant_file_path = os.path.join(work_path, "codon20.lst")
    # cmd_line = "Rscript " + kmeans_script_path + " " +  single_codon_file_path + " > " + cc_significant_file_path
    # os.system(cmd_line)

    # here we choose a  new way to infer the group 20140912



    # read the result of R process
    cc_significant_raw = RNAsnp.get_large_codon_group(cc_significant_file_path)
    cc_significant = [x for x in cc_significant_raw]
    print len(cc_significant)

    gu_model_file = "/home/zerodel/GitProjects/python-rna-structure/gu_model.mdl"
    rebuild_model_file = "/home/zerodel/GitProjects/python-rna-structure/rebuild_model.mdl"
    nested_model_file = "/home/zerodel/GitProjects/python-rna-structure/rna_structure.mdl"
    matrixh.degenerate(gu_model_file,
                       rebuild_model_file)

    matrixh.nest_export(rebuild_model_file,nested_model_file,cc_significant_file_path)


    #
    # CPD_site = get_whole_cpd(output_path)
    # CPDState = build_state_cpd(CPD_site, gc.codon_list_hypothesis, cc_significant)
    # print_state_cpd(CPDState,"60.def",codon_list_hypothesis)
    #
    # if not os.path.exists("CPD_site.db"):
    #     with open("CPD_site.db", "w") as writecpd:
    #         pickle.dump(CPD_site, writecpd)
    # else:
    #     with open("CPD_site.db", "r") as fcpd:
    #         CPD_site = pickle.load(fcpd)
    # print "----"
    # for key in CPDState.keys():
    #     if "psi" in CPDState[key]:
    #         if compareCodonPair(key)[1] == "T":
    # print key

if __name__ == "__main__":
   built_model()