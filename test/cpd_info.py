# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#
import os
import os.path
import pyRNAsnp.pyRNAsnp as RNAsnp
import pyHYPHY.global_constants as gc
import pyHYPHY.ModelHYPHY as pyMH

def test_check_codons():
    ccl_file = "/home/zerodel/Workspace/codon_lst"
    mp = pyMH.ModelHYPHY(lst_2_tuple_file_name=ccl_file)
    print str(mp.check_codons("AAC", "ACC"))



def cpd_reshape_export(cpd_path, ccl_file, output_file, codon_vector=gc.codon_list_hypothesis[:]):
    cpd_main = RNAsnp.get_whole_cpd(cpd_path)
    #codon_vector = gc.codon_list_hypothesis[:]
    mp = pyMH.ModelHYPHY(lst_2_tuple_file_name=ccl_file)
    with open(output_file, "w") as writer:
        writer.write("from\tto\tMutationCount\tMutationClass\n")
        available_cc = cpd_main.keys()
        for codon_origin in codon_vector:

            for codon_target in codon_vector:
                # build a line for lst file
                num_mutation = 0
                class_mutation = 0
                if (codon_origin, codon_target) in available_cc:
                    num_mutation = len(cpd_main[codon_origin, codon_target])
                    codon_check_response = mp.check_codons(codon_origin, codon_target)

                    if type(codon_check_response) == tuple:  # two codon only differ one site
                        is_transversion, is_sysnonmous, is_larger, nt_diff = codon_check_response
                        class_mutation += 1
                        if not is_sysnonmous:
                            class_mutation += 1
                        if is_larger:
                            class_mutation += 2

                writer.write("%d\t%d\t%d\t%d\n" % (codon_vector.index(codon_origin) + 1,
                                                   codon_vector.index(codon_target) + 1,
                                                   num_mutation,
                                                   class_mutation))

if __name__ == "__main__":
    # ccl_file = "/home/zerodel/Workspace/Yeast/codon_all.lst"
    # output_file = "/home/zerodel/GitProjects/python-rna-structure/data/YeastCPDInfo.list"
    # cpdpath = "/home/zerodel/Workspace/Yeast/yeast_cpd"
    # cpd_reshape_export(cpdpath, ccl_file, output_file)
    print "yeast ok "
    ccl_file_ecoli = "/home/zerodel/Workspace/codon_all.lst"
    output_file_ecoli = "/home/zerodel/GitProjects/python-rna-structure/data/EcoliCPDInfo.list"
    cpdpath_ecoli = "/home/zerodel/Workspace/cpd_store_site"
    cpd_reshape_export(cpdpath_ecoli, ccl_file_ecoli, output_file_ecoli)
    print "ecoli , ok"