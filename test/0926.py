# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme: move into yeast ....
#

import os
import os.path
import pyRNAsnp.pyRNAsnp as RNAsnp
import pyHYPHY.global_constants as gc

import pyHYPHY.MatrixHYPHY as matrixh
import pyHYPHY.BfHYPHY as BfH
import pyHYPHY.DataHYPHY as DH

tree_of_yeast = "((Smik:0.06505810202920048191,(Skud:0.06824775484329334563,(Sbay:0.06449537146489706108,(Sklu:0.44406379138184037814,Scas:0.34852255969430995242):0.27361297838594561549):0.02621412500451008112):0.02467797569334097968):0.02448288696591136016,Spar:0.03007100476252094409,Scer:0.04122423044401372916):0.0;"

list_yeast_species = ["Sbay", "Smik", "Skud", "Sklu", "Scas", "Spar", "Scer"]
file_yeast_aln_path = "/home/zerodel/Workspace/Yeast/yeast_aln/yeast_7_aln"
# input_file_folder = "/home/zerodel/Workspace/0927"
parent_path = "/home/zerodel/GitProject"
work_path = "/home/zerodel/Workspace/Yeast"
source_path = "/home/zerodel/Workspace/Yeast/yeast_rnasnp"
output_path = "/home/zerodel/Workspace/Yeast/yeast_cpd"
gu_model_file = "/home/zerodel/GitProjects/python-rna-structure/gu_model.mdl"
rebuild_model_file = "/home/zerodel/Workspace/Yeast/rebuild_model.mdl"
nested_model_file = "/home/zerodel/Workspace/Yeast/rna_structure_yeast_full_length.mdl"


def make_model():
    single_codon_file_path = os.path.join(work_path, "codon_lst")
    cc_significant_file_path = os.path.join(work_path, "codon_all.lst")

    #step 1 , set up  .cpd files

    rnasnp_files, snppathname = RNAsnp.snp_dir_list(source_path)
    RNAsnp.snp_dir_traversal(rnasnp_files, snppathname, output_path, True, True)
    # step 2 , print each codon lst file
    for codon in gc.codon_list_hypothesis:
        codon_with_direction = codon + "_"
        print "--->", codon_with_direction
        RNAsnp.cpd_dir_traversal(codon_with_direction, output_path, single_codon_file_path)

    # RNAsnp.group_divide(single_codon_file_path, cc_significant_file_path)

    matrixh.degenerate(gu_model_file,
                       rebuild_model_file)

    matrixh.nest_export(rebuild_model_file, nested_model_file, cc_significant_file_path)




def test_gene_species_name():
    aln_files = [single_file for single_file in os.listdir(file_yeast_aln_path) if ".aln" == single_file[-4:]]
    num_not_match = 0
    genes_for_7 = []
    print "there are ", str(len(aln_files)), "genes here"
    for file_aln_1 in aln_files:
        full_path = os.path.join(file_yeast_aln_path, file_aln_1)
        with open(full_path, "r") as reader:
            contents = reader.readlines()
            species = [line.split()[0].strip() for line in contents]
        is_7 = True
        for single_species in species:
            if not single_species in list_yeast_species:
                print "some gene not match in ", file_aln_1
                is_7 = False

        if is_7:
            print file_aln_1, "has 7 !"
            genes_for_7.append(file_aln_1.split(".")[0])
        else:
            num_not_match += 1

    print "In summary , ", str(num_not_match), "genes has less than 7 species"


def make_bf():
    os.chdir("/home/zerodel/Workspace/Yeast/0927/main_full_length")
    model_gtr = "rebuild_model.mdl"
    model_nest = "rna_structure_yeast_full_length.mdl"

    bf_maker_nest = BfH.HYPHYBatchFile(species_name="yeast",
                                      model_file=model_nest,
                                      bf_template_file="template_global")
    bf_maker_gtr = BfH.HYPHYBatchFile(species_name="yeast",
                                      model_file=model_gtr,
                                      bf_template_file="template_global")

    jobids = [file1.split(".")[0] for file1 in os.listdir(".") if ".aln" == file1[-4:]]

    for job_id in jobids:
        aln_file_path_full = os.path.join(file_yeast_aln_path, job_id + ".aln")
        input_file = os.path.join(os.curdir, job_id + ".input")

        if not os.path.exists(input_file):
            DH.aln2input(aln_file_path_full, input_file)

        bfgtr = job_id + "gtr.bf"
        bf_maker_gtr.write_batch_file(dot_input=job_id + ".input",
                                         dot_aln=job_id + ".aln",
                                         hyphy_result_file=job_id + "gtr.result",
                                         hyphy_batch_file=bfgtr)

        bfnest = job_id + "nest.bf"
        bf_maker_nest.write_batch_file(dot_input=job_id + ".input",
                                         dot_aln=job_id + ".aln",
                                         hyphy_result_file=job_id + "nest.result",
                                         hyphy_batch_file=bfnest)


def test_bf():
    os.chdir("/home/zerodel/Workspace/Yeast/0927/main_full_length")
    bf_files = [single_file for single_file in os.listdir(os.curdir) if ".bf" == single_file[-3:]]
    model_gtr = "rebuild_model.mdl"
    model_nest = "rna_structure_yeast_full_length.mdl"
    num_checked = 0
    for bf_file in bf_files:

        print str(num_checked)
        with open(bf_file, "r") as reader:
            content = reader.read()
        if "gtr" in bf_file:
            print "gtr bf"
            if model_nest in content:
                print "nest model in gtr bf", "bf_file"

        if "nest" in bf_file:
            print "nest bf"
            if model_gtr in content:
                print "gtr model in nest bf", "bf_file"
        num_checked += 1

if __name__ == "__main__":
    test_bf()
    print "all over"