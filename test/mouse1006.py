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
import pyHYPHY.DownStreamHYPHY as dsh
import pyHYPHY.MatrixHYPHY as matrixh
import pyHYPHY.DataHYPHY as DH
import pyHYPHY.BfHYPHY as BfH

# (((macaque:0.02602176453689536428,(chimp:0.00432283016370958121,human:0.00384313582540942496):0.00915641286880421179):0.03904538064229965549,(dog:0.06293344444628010126,(cow:0.07356214799052007702,horse:0.04838456774727731280):0.00508183725708111558):0.01752218773337325258):0.09364969675060515197,rat:0.03601197114550155204,mouse:0.02852390293940118213):0.0;

genes = "(((macaque,(chimp,human)),(dog,(cow,horse))),rat,mouse):0.0;"
def test_aln():
    aln_path = "/home/zerodel/Workspace/mouse/mouse_8_species"
    given_genes = "macaque,chimp,human,dog,cow,horse,rat,mouse".split(",")
    dsh.test_gene_species_match(aln_path, given_genes)


def mouse_prepare():
    parent_path = "/home/zerodel/GitProject"
    work_path = "/home/zerodel/Workspace/mouse"
    source_path = "/home/zerodel/Workspace/mouse/mouse_rnasnp"
    output_path = "/home/zerodel/Workspace/mouse/cpd"
    single_codon_file_path = os.path.join(work_path, "codon_codon_lst")
    cc_significant_file_path = os.path.join(work_path, "codon_all.lst")
    rnasnp_files, snppathname = RNAsnp.snp_dir_list(source_path)

    # # step 1 , set up  .cpd files
    # RNAsnp.snp_dir_traversal(rnasnp_files, snppathname, output_path, True)
    # # # step 2 , print each codon lst file
    # for codon in gc.codon_list_hypothesis:
    #     codon_with_direction = codon + "_"
    #     print "--->", codon_with_direction
    #     RNAsnp.cpd_dir_traversal(codon_with_direction, output_path, single_codon_file_path)

    RNAsnp.group_divide(single_codon_file_path, cc_significant_file_path)

    gu_model_file = "/home/zerodel/GitProjects/python-rna-structure/gu_model.mdl"
    rebuild_model_file = "/home/zerodel/GitProjects/python-rna-structure/rebuild_model.mdl"
    nested_model_file = os.path.join(work_path, "rna_structure_full_mouse.mdl")

    # matrixh.degenerate(gu_model_file,
    #                    rebuild_model_file)

    matrixh.nest_export(rebuild_model_file, nested_model_file, cc_significant_file_path)

def ccl_info(ccl_file):
    with open(ccl_file, "r") as ccl_fetch:
        contents = ccl_fetch.readlines()
    print ccl_file, '------\n'

    for line in contents:
        l_c = line.split()
        gene_id = l_c.pop(0)
        snp_num = len(l_c)

        print gene_id, "\t", str(snp_num)

def batch_for_batch():
    os.chdir("/home/zerodel/Workspace/mouse/FullLengthMouse")
    model_gtr = "rebuild_model.mdl"
    model_nest = "rna_structure_full_mouse.mdl"
    aln_path = "/home/zerodel/Workspace/mouse/mouse_8_species"

    bf_maker_nest = BfH.HYPHYBatchFile(species_name="mouse",
                                      model_file=model_nest,
                                      bf_template_file="template_global")
    bf_maker_gtr = BfH.HYPHYBatchFile(species_name="mouse",
                                      model_file=model_gtr,
                                      bf_template_file="template_global")

    jobids = [file1.split(".")[0] for file1 in os.listdir(aln_path) if ".aln" == file1[-4:]]

    for job_id in jobids:
        aln_file_path_full = os.path.join(aln_path, job_id + ".aln")
        input_file = os.path.join(os.curdir, job_id + ".input")

        if not os.path.exists(input_file):
            DH.aln2input(aln_file_path_full, input_file)

        bfgtr = job_id + "gtr.bf"
        bf_maker_gtr.write_batch_file(dot_input=job_id + ".input",
                                         dot_aln=aln_file_path_full,
                                         hyphy_result_file=job_id + "gtr.result",
                                         hyphy_batch_file=bfgtr)

        bfnest = job_id + "nest.bf"
        bf_maker_nest.write_batch_file(dot_input=job_id + ".input",
                                         dot_aln=aln_file_path_full,
                                         hyphy_result_file=job_id + "nest.result",
                                         hyphy_batch_file=bfnest)


if __name__ == "__main__":
    batch_for_batch()
