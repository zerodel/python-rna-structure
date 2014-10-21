# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#
__author__ = 'Zerodel'
import os
import shutil
import pyHYPHY.DataHYPHY as dh
import pyHYPHY.BfHYPHY as BfH

import pyHYPHY.MatrixHYPHY as matrixh
import pyRNAsnp.pyRNAsnp as RNAsnp
import pyHYPHY.global_constants as gc

def make_model():
    work_path = "d:/Workspace/Ecoli"
    source_path = "d:/Workspace/Ecoli/ecoli_snp"
    output_path = "d:/Workspace/Ecoli/cpd_full"
    single_codon_file_path = os.path.join(work_path, "codon_lst")
    cc_significant_file_path = os.path.join(work_path, "codon_all.lst")
        #step 1 , set up  .cpd files

    gu_model_file = "d:/Home/GitProject/python-rna-structure/gu_model.mdl"
    rebuild_model_file = "d:/Workspace/Ecoli/rebuild_model_full_length.mdl"
    nested_model_file = "d:/Workspace/Ecoli/rna_full_length_structure.mdl"
    rnasnp_files, snppathname = RNAsnp.snp_dir_list(source_path)
    RNAsnp.snp_dir_traversal(rnasnp_files, snppathname, output_path, True, True)
    ###step 2 , print each codon lst file
    for codon in gc.codon_list_hypothesis:
        codon_with_direction = codon + "_"
        print "--->", codon_with_direction
        RNAsnp.cpd_dir_traversal(codon_with_direction, output_path, single_codon_file_path)

    RNAsnp.group_divide(single_codon_file_path, cc_significant_file_path)

    matrixh.degenerate(gu_model_file,
                       rebuild_model_file)

    matrixh.nest_export(rebuild_model_file, nested_model_file, cc_significant_file_path)

def check_aln_full_species(folder_alns):
    all_alns = [single_file for single_file in os.listdir(folder_alns) if single_file[-4:] == ".aln"]
    aln_gene_num = []
    for single_file in all_alns:
        full_aln_path = os.path.join(folder_alns, single_file)
        genes, gene_lengths = dh.aln_info(full_aln_path)
        aln_gene_num.append(len(genes))

    max_num = max(aln_gene_num)
    aln_full_gene = []
    for indexI, aln in enumerate(all_alns):
        if aln_gene_num[indexI] == max_num:
            aln_full_gene.append(all_alns[indexI])

    return aln_full_gene


def main():
    # main part
    aln_files_folder = "d:/Workspace/Ecoli/ecoli_10_species"
    #aln_files_folder  = "d:/Workspace/Ecoli/test"
    target_path = "d:/Workspace/Ecoli/NoGapP1"

    aln_files = [single_file for single_file in os.listdir(aln_files_folder) if ".aln" == os.path.splitext(single_file)[-1]]

    model_gtr = "rebuild_model.mdl"
    model_nest = "rna_full_length_structure.mdl"

    job_ids = []
    seq_length = []
    for aln in aln_files:
        aln_path_full = os.path.join(aln_files_folder,aln)

        gene_in_aln, length_in_aln = dh.aln_info(aln_path_full)
        seq_length.append(length_in_aln[0])
        job_id = aln.split(".")[0]
        if 10 == len(gene_in_aln): #only those gene shared in 10 species
            job_ids.append(job_id)
            input_file_path_full = os.path.join(target_path, job_id + ".input")
            shutil.copyfile(aln_path_full, os.path.join(target_path,aln))
            # write input file
            dh.aln2inputNogap(aln_path_full, input_file_path_full)

            # make bf file for full length no gap gtr and empirical model

    os.chdir(target_path)
    bf_maker_nest = BfH.HYPHYBatchFile(species_name="ecoli",
                                      model_file=model_nest,
                                      bf_template_file="templateNoGap")
    bf_maker_gtr = BfH.HYPHYBatchFile(species_name="ecoli",
                                      model_file=model_gtr,
                                      bf_template_file="templateNoGap")

    for job_id in job_ids:
        # aln_path_full = os.path.join(aln_files_folder,aln)
        # input_file_path_full = os.path.join(target_path, job_id + ".input")
        # shutil.copyfile(aln_path_full, os.path.join(target_path,aln))
        # # write input file
        # dh.aln2inputNogap(aln_path_full, input_file_path_full)

        # make bf file for full length no gap gtr and empirical model

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

    bf_maker_gtr.set_partition(0,51)
    bf_maker_nest.set_partition(0,51)

    for job_id in job_ids:
        bfgtr = job_id + "gtrp1.bf"
        bf_maker_gtr.write_batch_file(dot_input=job_id + ".input",
                                     dot_aln=job_id + ".aln",
                                     hyphy_result_file=job_id + "gtrp1.result",
                                     hyphy_batch_file=bfgtr)

        bfnest = job_id + "nestp1.bf"
        bf_maker_nest.write_batch_file(dot_input=job_id + ".input",
                                     dot_aln=job_id + ".aln",
                                     hyphy_result_file=job_id + "nestp1.result",
                                     hyphy_batch_file=bfnest)


    for indexI, job_id in enumerate(job_ids):
        bf_maker_gtr.set_partition(51, seq_length[indexI])
        bf_maker_nest.set_partition(51, seq_length[indexI])

        bfgtr = job_id + "gtrp2.bf"
        bf_maker_gtr.write_batch_file(dot_input=job_id + ".input",
                                     dot_aln=job_id + ".aln",
                                     hyphy_result_file=job_id + "gtrp2.result",
                                     hyphy_batch_file=bfgtr)

        bfnest = job_id + "nestp2.bf"
        bf_maker_nest.write_batch_file(dot_input=job_id + ".input",
                                     dot_aln=job_id + ".aln",
                                     hyphy_result_file=job_id + "nestp2.result",
                                     hyphy_batch_file=bfnest)



if __name__ == "__main__":
    main()