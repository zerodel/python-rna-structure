# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#
__author__ = 'Zerodel'
import os
import os.path
import pyHYPHY.BfHYPHY as bf
import pyHYPHY.global_constants as gc
def test_bf():
    bf1 = bf.HYPHYBatchFile(species_name="ecoli",model_file="../rebuild_model.mdl",bf_template_file="../pyHYPHY/templateNoGap")
    bf1.use_given_tree = True
    bf1.tree_definition_external = "(g1,g2)"

    bf1.write_batch_file("./tmp2.input", "./tmp.aln","xxx.bf","xxx.result")



if __name__ == "__main__":
    test_bf()