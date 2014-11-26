# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#
__author__ = 'Zerodel'

import os
import os.path

def AddSysPath(new_path):
    """ AddSysPath(new_path): adds a "directory" to Python's sys.path
    Does not add the directory if it does not exist or if it's already on
    sys.path. Returns 1 if OK, -1 if new_path does not exist, 0 if it was
    already on sys.path.
    """
    import sys, os
    # Avoid adding nonexistent paths
    if not os.path.exists(new_path): return -1
    # Standardize the path.  Windows is case-insensitive, so lowercase
    # for definiteness if we are on Windows.
    new_path = os.path.abspath(new_path)
    if sys.platform == 'win32':
        new_path = new_path.lower( )
    # Check against all currently available paths
    for x in sys.path:
        x = os.path.abspath(x)
        if sys.platform == 'win32':
            x = x.lower( )
        if new_path in (x, x + os.sep):
            return 0
    sys.path.append(new_path)
    # if you want the new_path to take precedence over existing
    # directories already in sys.path, instead of appending, use:
    # sys.path.insert(0, new_path)
    return 1

AddSysPath("..")


import pyHYPHY.DataHYPHY as DH
import pyHYPHY.BfHYPHY as BfH
import pyHYPHY.SpeciesSpecificTree as SS

class AlignmentNotSame(Exception):
    pass

class DimNotSame(Exception):
    pass

class SequenceTooShort(Exception):
    pass


# first , check  the gene  order in all alignment files
def check_gene_order_in_alignments(path_to_alignments):
    """
    判断path_to_alignment 里面所有的.aln 文件里物种名是否都是一样的.

    :param path_to_alignments:
    :return:
    """
    # get file-names
    import os
    import os.path

    current_dir = os.path.abspath(os.curdir)
    os.chdir(path_to_alignments)
    aln_files = [single_file for single_file in os.listdir(path_to_alignments)
                 if ".aln" == os.path.splitext(single_file)[-1]]

    try:
        import pyHYPHY.DataHYPHY as dh

    except ImportError:
        print "error in importing pyHYPHY module"
        return

    aln_genes_calibration = dh.aln_info(aln_files[0])[0]
    gene_un_match = 0
    for aln_entry in aln_files:
        if not aln_genes_calibration == dh.aln_info(aln_entry)[0]:
            print aln_entry, ":---- ", str(dh.aln_info(aln_entry)[0]), "\n"
            gene_un_match += 1

    print "whole number of unmatch file is %d \n and pattern is :\n %s \n" % (gene_un_match, aln_genes_calibration)

    os.chdir(current_dir)

    if 0 == gene_un_match:
        return aln_genes_calibration
    else:
        raise AlignmentNotSame


def paste_matrix(matrix_main, matrix_tail):
    """
    链接两个字符串数组.

    两个字符串数组元素个数不一样, 将出现DimNotSame异常.

    :param matrix_main: 被连接的字符串数组
    :param matrix_tail: 添加在matrix_main 后面的字符串数组.
    :return: 拼接好的字符串数组
    """
    if not len(matrix_main) == len(matrix_tail):
        raise DimNotSame
    tmp = []
    for indexI, line in enumerate(matrix_main):
        tmp.append("".join([line, matrix_tail[indexI]]))

    return tmp


def extract_TIR_single_file(aln_file, length_TIR, start_point=0):
    """
    读取aln_file,把序列start_point 开始的 length_TIR 部分提取出来.

    :param aln_file:
    :param length_TIR:
    :param start_point:
    :return: 一个字符串数组.
    """
    with open(aln_file, "r") as reader:
        sequences = [line.split()[-1] for line in reader.readlines()]

    if len(sequences[0]) < start_point+length_TIR:
        raise SequenceTooShort

    matrix_TIR_raw = ["".join(line[start_point:start_point+length_TIR]) for line in sequences]

    return matrix_TIR_raw


<<<<<<< HEAD
def gene_conjunction(path_aln_file, input_file_for_hyphy, length_of_TIR, start_point=0):
=======
def gene_conjunction(path_aln_file, input_file_for_hyphy, length_of_TIR, start_point=0, remove_gap=True):
>>>>>>> f0acd743d1106c96b88083b2df2cb3526b388aec
    # combination of all translation initiation region

    import os
    import os.path

    current_dir = os.path.abspath(os.curdir)
    os.chdir(path_aln_file)
    species_list = check_gene_order_in_alignments(path_aln_file)
    inputfile_header = [">%s" % species_name for species_name in species_list]
    matrix_sequence = ["" for species_name in species_list]

    aln_files = [single_file for single_file in os.listdir(path_aln_file)
                 if ".aln" == os.path.splitext(single_file)[-1]]

    for single_aln_file in aln_files:
        # single file operation
        # rejection : 1. length not enough 2 two many  gaps
        try:
<<<<<<< HEAD
            matrix_sequence = paste_matrix(matrix_sequence,DH.remove_gaps_matrix(extract_TIR_single_file(single_aln_file, length_of_TIR, start_point)))
=======

            if remove_gap:
                gene_seq_matrix_addition = DH.remove_gaps_matrix(extract_TIR_single_file(single_aln_file, length_of_TIR, start_point))
            else:
                gene_seq_matrix_addition = extract_TIR_single_file(single_aln_file, length_of_TIR, start_point)

            matrix_sequence = paste_matrix(matrix_sequence, gene_seq_matrix_addition)
>>>>>>> f0acd743d1106c96b88083b2df2cb3526b388aec

        except SequenceTooShort:
            print "Too short in %s" % single_aln_file
            continue
        except DimNotSame:
            print "Error of Dimisions %s" % single_aln_file
            continue
        else:
            pass

    with open(input_file_for_hyphy, "w") as writerhere:
        print "input file : %s in %s" % (input_file_for_hyphy, os.path.abspath(os.curdir))
        for indexI, gene_title in enumerate(inputfile_header):
            writerhere.write(gene_title + "\n")
            writerhere.write(matrix_sequence[indexI] + "\n")

    os.chdir(current_dir)

def make_bf_TIR():
    workpath = "d:\Workspace\Ecoli\combinedP1"
    os.chdir(workpath)
    model_gtr = "rebuild_model.mdl"
    model_nest = "rna_full_length_structure.mdl"

    bf_maker_nest = BfH.HYPHYBatchFile(species_name="ecoli",
                                       model_file=model_nest,
                                       bf_template_file="templateNoGap")
    bf_maker_gtr = BfH.HYPHYBatchFile(species_name="ecoli",
                                      model_file=model_gtr,
                                      bf_template_file="templateNoGap")

    bf_maker_gtr.set_tree_from_outside(SS.ecoli())
    bf_maker_nest.set_tree_from_outside(SS.ecoli())

    for window_width_offset in range(15):
        length_window = 30 + 5*window_width_offset
        aln_files = "d:\Workspace\Ecoli\ecoli_10_species"

        extracted_input_file = workpath + "\TIR" + str(length_window) + ".input"
        gene_conjunction(aln_files, extracted_input_file, length_window)

        bf_maker_gtr.write_batch_file(dot_input="TIR%d.input" % length_window,
                                      dot_aln="",
                                      hyphy_batch_file="TIR%dgtr.bf" % length_window,
                                      hyphy_result_file="TIR%dgtr.result" % length_window)

        bf_maker_nest.write_batch_file(dot_input="TIR%d.input" % length_window,
                                       dot_aln="",
                                       hyphy_batch_file="TIR%dnest.bf" % length_window,
                                       hyphy_result_file="TIR%dnest.result" % length_window)


def make_bf_p2():
    workpath = "d:\Workspace\Ecoli\P2"
    os.chdir(workpath)
    model_gtr = "rebuild_model.mdl"
    model_nest = "rna_full_length_structure.mdl"

    bf_maker_gtr = BfH.HYPHYBatchFile(species_name="ecoli",
                                      model_file=model_gtr,
                                      bf_template_file="templateNoGap")

    bf_maker_nest = BfH.HYPHYBatchFile(species_name="ecoli",
                                       model_file=model_nest,
                                       bf_template_file="templateNoGap")

    bf_maker_gtr.set_tree_from_outside(SS.ecoli())
    bf_maker_nest.set_tree_from_outside(SS.ecoli())

    aln_file_folder = "d:\Workspace\Ecoli\ecoli_10_species"

    aln_files = [file_single for file_single in os.listdir(aln_file_folder) if ".aln" == os.path.splitext(file_single)[-1]]

    for single_aln in aln_files:
        aln_full_path = os.path.join(aln_file_folder, single_aln)
        genes, lengths = DH.aln_info(aln_full_path)
        gene_full_length = lengths[0]
        jobid = single_aln.split(".")[0]
        input_file_name = "%s.input" % jobid
        DH.aln2input(aln_full_path, input_file_name)

        if gene_full_length < 52:
            print "%s too short ---" % single_aln
            continue

        bf_maker_gtr.set_partition(51, gene_full_length)
        bf_maker_nest.set_partition(51, gene_full_length)

        bf_maker_gtr.write_batch_file(dot_input=input_file_name,
                                      dot_aln="",
                                      hyphy_batch_file="%sp2gtr.bf" % jobid,
                                      hyphy_result_file="%sp2gtr.result" % jobid)


        bf_maker_nest.write_batch_file(dot_input=input_file_name,
                                       dot_aln="",
                                       hyphy_batch_file="%sp2nest.bf" % jobid,
                                       hyphy_result_file="%sp2nest.result" % jobid)


def sliding_window():
    """perform a sliding window analysis over sequence alignment
"""
    workpath = "/Users/zerodel/WorkSpace/test/test_slidingWindow"

    # important here, has changed path.
    os.chdir(workpath)
    model_gtr = "rebuild_model.mdl"
    model_nest = "rna_full_length_structure.mdl"
    aln_files = "/Users/zerodel/WorkSpace/ecoli_aln"

    bf_maker_nest = BfH.HYPHYBatchFile(species_name="ecoli",
                                       model_file=model_nest,
                                       bf_template_file="templateNoGap")
    bf_maker_gtr = BfH.HYPHYBatchFile(species_name="ecoli",
                                      model_file=model_gtr,
                                      bf_template_file="templateNoGap")

    bf_maker_gtr.set_tree_from_outside(SS.ecoli())
    bf_maker_nest.set_tree_from_outside(SS.ecoli())

    step_num = 21
    step_width_step = 15
    window_width_min = 30
    site_start_point = 0
    site_shift_step = 15

    window_axis = [0,2]
    start_axis = range(step_num)

    for window_width_offset in window_axis:
        for window_start_offset in start_axis:

            start_site = site_shift_step*window_start_offset + site_start_point
            length_window = step_width_step*window_width_offset + window_width_min

            # extracted_input_file = workpath + "\TIR" + str(length_window) + ".input"
            job_description = "w%ds%d" % (length_window, start_site)
            extracted_input_file = "TIR%s.input" % job_description
            extracted_input_no_gap = "TIR%sN.input" % job_description

            gene_conjunction(aln_files, os.path.join(workpath, extracted_input_file), length_window, start_site)
            gene_conjunction(aln_files, os.path.join(workpath, extracted_input_no_gap), length_window, start_site, False)


            bf_maker_gtr.write_batch_file(dot_input=extracted_input_file,
                                          dot_aln="",
                                          hyphy_batch_file="TIR%sgtr.bf" % job_description,
                                          hyphy_result_file="TIR%sgtr.result" % job_description)

            bf_maker_nest.write_batch_file(dot_input=extracted_input_file,
                                           dot_aln="",
                                           hyphy_batch_file="TIR%snest.bf" % job_description,
                                           hyphy_result_file="TIR%snest.result" % job_description)



                # no gap
            bf_maker_gtr.write_batch_file(dot_input=extracted_input_no_gap,
                                          dot_aln="",
                                          hyphy_batch_file="TIR%sgtrN.bf" % job_description,
                                          hyphy_result_file="TIR%sgtrN.result" % job_description)

            bf_maker_nest.write_batch_file(dot_input=extracted_input_no_gap,
                                           dot_aln="",
                                           hyphy_batch_file="TIR%snestN.bf" % job_description,
                                           hyphy_result_file="TIR%snestN.result" % job_description)

def sliding_window():
    workpath = "d:\Workspace\Ecoli\slidingWindow"
    os.chdir(workpath)
    model_gtr = "rebuild_model.mdl"
    model_nest = "rna_full_length_structure.mdl"

    bf_maker_nest = BfH.HYPHYBatchFile(species_name="ecoli",
                                      model_file=model_nest,
                                      bf_template_file="templateNoGap")
    bf_maker_gtr = BfH.HYPHYBatchFile(species_name="ecoli",
                                      model_file=model_gtr,
                                      bf_template_file="templateNoGap")

    bf_maker_gtr.set_tree_from_outside(SS.ecoli())
    bf_maker_nest.set_tree_from_outside(SS.ecoli())

    step_num = 15
    step_width_step = 5
    window_width_min = 30
    site_start_point = 0
    site_shift_step = 5
    for window_start_offset in range(step_num):
        for window_width_offset in range(step_num):

            start_site = site_shift_step*window_start_offset + site_start_point
            length_window = step_width_step*window_width_offset + window_width_min
            aln_files = "d:\Workspace\Ecoli\ecoli_10_species"

            # extracted_input_file = workpath + "\TIR" + str(length_window) + ".input"
            job_description = "w%ds%d" % (length_window, start_site)
            extracted_input_file = "%s\\TIR%s.input" % (workpath, job_description)
            print extracted_input_file

            gene_conjunction(aln_files, extracted_input_file, length_window, start_site)

            bf_maker_gtr.write_batch_file(dot_input="TIR%s.input" % job_description,
                                          dot_aln="",
                                          hyphy_batch_file="TIR%sgtr.bf" % job_description,
                                          hyphy_result_file="TIR%sgtr.result" % job_description)

            bf_maker_nest.write_batch_file(dot_input="TIR%s.input" % job_description,
                                           dot_aln="",
                                          hyphy_batch_file="TIR%snest.bf" % job_description,
                                          hyphy_result_file="TIR%snest.result" % job_description)


def make_bf_p2():
    workpath = "d:\Workspace\Ecoli\P2"
    os.chdir(workpath)
    model_gtr = "rebuild_model.mdl"
    model_nest = "rna_full_length_structure.mdl"

    bf_maker_gtr = BfH.HYPHYBatchFile(species_name="ecoli",
                                      model_file=model_gtr,
                                      bf_template_file="templateNoGap")

    bf_maker_nest = BfH.HYPHYBatchFile(species_name="ecoli",
                                      model_file=model_nest,
                                      bf_template_file="templateNoGap")

    bf_maker_gtr.set_tree_from_outside(SS.ecoli())
    bf_maker_nest.set_tree_from_outside(SS.ecoli())

    aln_file_folder = "d:\Workspace\Ecoli\ecoli_10_species"

    aln_files = [file_single for file_single in os.listdir(aln_file_folder) if ".aln" == os.path.splitext(file_single)[-1]]


    for single_aln in aln_files:
        aln_full_path = os.path.join(aln_file_folder, single_aln)
        genes, lengths = DH.aln_info(aln_full_path)
        gene_full_length = lengths[0]
        jobid = single_aln.split(".")[0]
        input_file_name = "%s.input" % jobid
        DH.aln2input(aln_full_path, input_file_name)

        if gene_full_length < 52:
            print "%s too short ---" % single_aln
            continue

        bf_maker_gtr.set_partition(51, gene_full_length)
        bf_maker_nest.set_partition(51, gene_full_length)

        bf_maker_gtr.write_batch_file(dot_input=input_file_name,
                                      dot_aln="",
                                      hyphy_batch_file="%sp2gtr.bf" % jobid,
                                      hyphy_result_file="%sp2gtr.result" % jobid)


        bf_maker_nest.write_batch_file(dot_input=input_file_name,
                                      dot_aln="",
                                      hyphy_batch_file="%sp2nest.bf" % jobid,
                                      hyphy_result_file="%sp2nest.result" % jobid)



if __name__ == "__main__":
<<<<<<< HEAD
    print "hello world"
=======
    print SS.ecoli()
>>>>>>> f0acd743d1106c96b88083b2df2cb3526b388aec
