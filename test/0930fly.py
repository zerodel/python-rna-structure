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

# (((dyak:0.02222387244577254603,dere:0.02474897373042929741):0.00750822541146930773,(((dpse:0.00322057838860367309,dper:0.00795217445289250757):0.12008827375416616934,((dgri:0.11476940268550095414,(dmoj:0.10284627994071587898,dvir:0.06466471721674015016):0.02810781121689204590):0.11666583477695514903,dwil:0.17403072889350790220):0.05051801556694598960):0.05346497847295386685,dana:0.10682155767321652173):0.08080725121552659318):0.01626751927252095831,(dsec:0.00888712911143530772,dsim:0.00840237751560124480):0.00824132139066967133,dmel:0.01370395896674640214):0.0;
#dyak,dere, dpse,dper,dgri,dmoj,dvir,dwil,dana,dsec,dsim,dmel


def fly_prepare():
    parent_path = "/home/zerodel/GitProject"
    work_path = "/home/zerodel/Workspace/Fly"
    source_path = "/home/zerodel/Workspace/Fly/flyrnasnp"
    output_path = "/home/zerodel/Workspace/Fly/cpd"
    single_codon_file_path = os.path.join(work_path, "codon_codon_lst")
    cc_significant_file_path = os.path.join(work_path, "codon_all.lst")
    rnasnp_files, snppathname = RNAsnp.snp_dir_list(source_path)

    # # step 1 , set up  .cpd files
    # RNAsnp.snp_dir_traversal(rnasnp_files, snppathname, output_path, True)
    # # step 2 , print each codon lst file
    # for codon in gc.codon_list_hypothesis:
    #     codon_with_direction = codon + "_"
    #     print "--->", codon_with_direction
    #     RNAsnp.cpd_dir_traversal(codon_with_direction, output_path, single_codon_file_path)

    # RNAsnp.group_divide(single_codon_file_path, cc_significant_file_path)

    gu_model_file = "/home/zerodel/GitProjects/python-rna-structure/gu_model.mdl"
    rebuild_model_file = "/home/zerodel/GitProjects/python-rna-structure/rebuild_model.mdl"
    nested_model_file = os.path.join(work_path, "rna_structure_full_fly.mdl")

    # matrixh.degenerate(gu_model_file,
    #                    rebuild_model_file)

    matrixh.nest_export(rebuild_model_file, nested_model_file, cc_significant_file_path)


def test_aln():
    aln_path = "/home/zerodel/Workspace/Fly/align_fly/fly_12_species"
    given_genes = "dyak,dere,dpse,dper,dgri,dmoj,dvir,dwil,dana,dsec,dsim,dmel".split(",")
    dsh.test_gene_species_match(aln_path, given_genes)


def make_bf():
    os.chdir("/home/zerodel/Workspace/Fly/full_length_0930")
    model_gtr = "rebuild_model.mdl"
    model_nest = "rna_structure_full_fly.mdl"
    aln_path = "/home/zerodel/Workspace/Fly/align_fly/fly_12_species"

    bf_maker_nest = BfH.HYPHYBatchFile(species_name="fly",
                                      model_file=model_nest,
                                      bf_template_file="template_global")
    bf_maker_gtr = BfH.HYPHYBatchFile(species_name="fly",
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
    # fly_prepare()
    # test_aln()
    make_bf()