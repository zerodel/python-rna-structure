# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#
__author__ = 'Zerodel'

import os
import os.path
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
    if not len(matrix_main) == len(matrix_tail):
        raise DimNotSame
    tmp = []
    for indexI, line in enumerate(matrix_main):
        tmp.append("".join([line, matrix_tail[indexI]]))

    return tmp



def extract_TIR_single_file(aln_file, length_TIR):
    with open(aln_file, "r") as reader:
        sequences = [line.split()[-1]  for line in reader.readlines()]

    if len(sequences[0]) < length_TIR:
        raise SequenceTooShort

    matrix_TIR_raw = ["".join(line[0:length_TIR]) for line in sequences]

    return matrix_TIR_raw


def gene_conjunction(path_aln_file, input_file_for_hyphy, length_of_TIR):
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
        # rejection : 1. length not enough 2 two many gaps
        try:
            matrix_sequence = paste_matrix(matrix_sequence,DH.remove_gaps_matrix(extract_TIR_single_file(single_aln_file, length_of_TIR)))

        except SequenceTooShort:
            print "Too short in %s" % single_aln_file
            continue
        except DimNotSame:
            print "Error of Dimisions %s" % single_aln_file
            continue

        else:
            pass


    with open(input_file_for_hyphy, "w") as writerhere:
        for indexI , gene_title in enumerate(inputfile_header):
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


if __name__ == "__main__":
    make_bf_TIR()