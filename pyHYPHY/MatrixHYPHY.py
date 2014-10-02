# !/usr/bin/env python
# -*- coding : utf-8 -*-
# author :
# Readme:
# changelog :nothing
#
__author__ = 'Zerodel_2'

import re
import pyHYPHY.global_constants as gc
import pyHYPHY.ModelHYPHY as pyMH

import logging
logging.basicConfig(filename="MatrixHYPHY.log", filemode="w", level=logging.DEBUG)


def degenerate(gtr_nested_file, gtr_mdl_file):
    """deduce the GTR model from Gu's 2010 model
    :return:
    """

    with open(name=gtr_nested_file, mode="r") as readerG:
        model_gu = [line for line in readerG if "{" in line or "}" in line]

    rebuilt_lines = []
    psi_removed = 0
    for line_num, line in enumerate(model_gu):
        entry_in_line = line.split(",")
        newline = []
        for entry_index, single_entry in enumerate(entry_in_line):
            if "*psi" in single_entry:
                # remove "*psi"
                new_entry = re.sub("\*psi", "", single_entry)
                if not new_entry == single_entry:
                    psi_removed += 1
            else:
                new_entry = single_entry

            newline.append(new_entry)

        line_rebuilt = ",".join(newline)  # rebuilt a line
        rebuilt_lines.append(line_rebuilt)
    model_rebuilt = "".join(rebuilt_lines)  # rebuilt a model
    print "%d psi removed" % psi_removed
    with open(name=gtr_mdl_file, mode="w") as model_w:
        model_w.write(model_rebuilt)


def nest_export(gtr_mdl, gtr_nested, ccl_file):
    """
    export a nested model, need a file contains the classic gtr model , and a file contains CCL file

    Args:
    gtr_mdl : file name of existing gtr model
    gtr_nested : file name for output nested model
    CCL_file : a file contains list of codon-codon pairs , which assigned the psi parameter

    Return:

    """
    with open(name=gtr_mdl, mode="r") as readerG:
        model_gu = [line for line in readerG if "{" in line or "}" in line]

    index_codon = gc.codon_list_hypothesis[:]

    psi_in_model = 0
    end_of_line = 0
    mp = pyMH.ModelHYPHY(lst_2_tuple_file_name=ccl_file)
    rebuilt_lines = []

    for line_num, line in enumerate(model_gu):
        entry_in_line = line.split(",")
        newline = []
        for entry_index, single_entry in enumerate(entry_in_line):
            new_entry = single_entry
            codon_check_response = mp.check_codons(codon_origin=index_codon[line_num],
                                                   codon_target=index_codon[entry_index])

            if type(codon_check_response) == tuple:  # two codon only differ one site
                is_transversion, is_sysnonmous, is_larger, nt_diff = codon_check_response
                if is_larger:

                    if "}" in single_entry:
                        end_of_line += 1
                        pos = single_entry.index("}")
                        new_entry = single_entry[0:pos] + "*psi" + single_entry[pos:-1] + "\n"
                    else:
                        new_entry = single_entry + "*psi"
                logging.debug("\n codon1:%s \t codon2:%s \t entry:%s \n" %
                              (index_codon[line_num], index_codon[entry_index], new_entry))

            newline.append(new_entry)
            if "*psi" in new_entry:
                psi_in_model += 1

        line_rebuilt = ",".join(newline)  # rebuilt a line
        rebuilt_lines.append(line_rebuilt)
    model_rebuilt = "".join(rebuilt_lines)  # rebuilt a model

    print end_of_line, " ends in model"
    with open(name=gtr_nested, mode="w") as model_w:
        model_w.write(model_rebuilt)

    return model_rebuilt

if __name__ == "__main__":
    pass