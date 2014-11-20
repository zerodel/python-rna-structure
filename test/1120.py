# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#

import sys
import os
import os.path

par_path = os.path.abspath(os.pardir)
if not par_path in sys.path:
    sys.path.append(par_path)

# assume all those files contain same number of species, for example : 10 
class InputDict(dict):
    def paste(self, aln_line):
        gene_id, gene_seq = aln_line.split()
        gene_id = gene_id.strip()
        self.setdefault(gene_id,[])
        gene_seq = self._cleaned(gene_seq)
        self[gene_id].append(gene_seq)
    
    def _cleaned(self, gene_seq):
        return self._sequence_trim(gene_seq)

    def _sequence_trim(self, gene_seq):
        return gene_seq[51:]
######################################################################################
### for test in emacs python mode
test_emacs = True
if __name__ == "__main__":
    test_emacs = False

if (test_emacs):
    import random
    import os
    import os.path
    import shutil

    lowpsi = [id_aln.strip() for id_aln in "b0440 b1823  b3313 b3321 b3736".split()]
    highpsi = [id_aln.strip() for id_aln in "b0096 b0171  b0471 b0884 b2235 b3065 b3185  b3302 b3783".split()]
    source_path = "/media/zerodel/新加卷2/Workspace/Ecoli/NoGapP1"
    lowpsi_path = "/media/zerodel/新加卷2/Workspace/Ecoli/test/lowpsi"
    highpsi_path = "/media/zerodel/新加卷2/Workspace/Ecoli/test/highpsi"


    for psi_id in lowpsi:
        aln_file_name = os.path.join(source_path, psi_id + ".aln")
        shutil.copy(aln_file_name, lowpsi_path + "/combined")
        # input_name = os.path.join(source_path, psi_id + ".input")
        # shutil.copy(input_name, lowpsi_path)
        # bf_name_nest = os.path.join(source_path, psi_id + "nest.bf")
        # shutil.copy(bf_name_nest, lowpsi_path)
        # bf_name_gtr = os.path.join(source_path, psi_id + "gtr.bf")
        # shutil.copy(bf_name_gtr, lowpsi_path)

    for psi_id in highpsi:
        aln_file_name = os.path.join(source_path, psi_id + ".aln")
        shutil.copy(aln_file_name, highpsi_path + "/combined")
        # input_name = os.path.join(source_path, psi_id + ".input")
        # shutil.copy(input_name, highpsi_path)
        # bf_name_nest = os.path.join(source_path, psi_id + "nest.bf")
        # shutil.copy(bf_name_nest, highpsi_path)
        # bf_name_gtr = os.path.join(source_path, psi_id + "gtr.bf")
        # shutil.copy(bf_name_gtr, highpsi_path)
# combination part
    seq_matrix = InputDict()
    workpath = "/media/zerodel/新加卷2/Workspace/Ecoli/test/highpsi"
    aln_path = workpath + "/combined"
    alnfiles = [os.path.join(aln_path, file_1) for file_1 in os.listdir(aln_path)
                if ".aln" == os.path.splitext(file_1)[-1]]

    file_picked = random.choice(alnfiles)
    print file_picked
    gene_order = []
    with open(file_picked) as reader1:
        for line in reader1:
            gene_order.append(line.split()[0].strip())

    print gene_order

    for file1 in alnfiles:
        with open(file1) as reader2:
            for line in reader2:
              seq_matrix.paste(line)  
    
    print [[len(seq) for seq in seq_matrix[gene_name]] for gene_name in gene_order]
    with open(os.path.join(aln_path, "Tmp.input"), "w") as writer:
        for gene_id in gene_order:
            writer.write(">%s\n" % gene_id)
            writer.write("%s\n" % "".join(seq_matrix[gene_id]))

If __name__ == "__main__":
    pass
    # here to add some self desciption.....
