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




def main1():
    cpd_main = RNAsnp.get_whole_cpd("/home/zerodel/Workspace/cpd_store_site")
    codon_vector = gc.codon_list_hypothesis
    with open("/home/zerodel/GitProjects/python-rna-structure/data/rnasnp.matrix", "w") as writer:
        for codon1 in codon_vector:
            # build a line for lst file
            line_single = []
            for codon2 in codon_vector:
                line_single.append(str(len(cpd_main[(codon1, codon2)])))

            writer.write("\t".join(line_single) + "\n")
    with open("/home/zerodel/GitProjects/python-rna-structure/data/codon_vector.txt", "w") as writer:
        writer.write("\n".join(gc.codon_list_hypothesis))

if __name__ == "__main__":
    parent_path = "/home/zerodel/GitProject"
    work_path = "/home/zerodel/Workspace"
    source_path = "/home/zerodel/Workspace/yeast_rnasnp"
    output_path = "/home/zerodel/Workspace/yeast_cpd"
    single_codon_file_path = os.path.join(work_path, "codon_lst")
    cc_significant_file_path = os.path.join(work_path, "codon_all.lst")
    rnasnp_files, snppathname = RNAsnp.snp_dir_list(source_path)

    # #step 1 , set up  .cpd files
    # RNAsnp.snp_dir_traversal(rnasnp_files, snppathname, output_path, True)
    cpd_main = RNAsnp.get_whole_cpd(output_path)
    codon_vector = gc.codon_list_64
    with open("/home/zerodel/GitProjects/python-rna-structure/data/YeastfullRnasnp.matrix", "w") as writer:
        available_cc = cpd_main.keys()
        for codon1 in codon_vector:
            # build a line for lst file
            line_single = []
            for codon2 in codon_vector:
                if (codon1, codon2) in available_cc:
                    line_single.append(str(len(cpd_main[(codon1, codon2)])))
                else:
                    line_single.append(str(0))

            writer.write("\t".join(line_single) + "\n")