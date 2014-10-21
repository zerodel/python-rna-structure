# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author :
# Readme: for all  . aln files in one folder , export a .input file and process four model
# changelog :

__author__ = 'Zerodel_2'

import os
import os.path
import pyHYPHY.BfHYPHY as bfHYPHY
import pyHYPHY.DataHYPHY as dataHYPHY


def aln_folder_traversal(folder_name):
    """    """
    built_in_gy94mdl = "/usr/lib/hyphy/TemplateBatchFiles/TemplateModels/GY94.mdl"

    gy94bf = bfHYPHY.HYPHYBatchFile(species_name="ecoli",
                                    model_file=built_in_gy94mdl,
                                    bf_template_file="partition.bf")

    own_model = "nest_gy.mdl"
    nested_model = bfHYPHY.HYPHYBatchFile(species_name="ecoli",
                                          model_file=own_model,
                                          bf_template_file="partition.bf")

    nt_nest_model = "nt_nest_gy.mdl"
    bf_nt_nest = bfHYPHY.HYPHYBatchFile(species_name="ecoli",
                                     model_file=nt_nest_model,
                                     bf_template_file="partition.bf")

    gu_model = "myCodonMatrix.def"
    gu_bf = bfHYPHY.HYPHYBatchFile(species_name="ecoli",
                                      model_file=gu_model,
                                      bf_template_file="synAlphaWPsiModelP.bf")

    pwd = os.path.abspath(folder_name)
    os.chdir(folder_name)
    aln_files = [file1
                 for file1 in os.listdir(pwd)
                 if "aln" == file1.split(".")[-1]]

    for index, aln_name in enumerate(aln_files):
        # write batch file for each aln
        gene_id = aln_name.split(".")[0]
        genes, gene_len = dataHYPHY.aln_info(aln_name)
        len_gene = max(gene_len)

        gy94bf1 = "%sgy94.bf" % gene_id
        input_filename = "%s.input" % gene_id

        # input file here
        dataHYPHY.aln2input(dot_aln_file=aln_name, hyphy_input_file=input_filename)

        # gy model , full length
        gy94bf.write_batch_file(dot_aln=aln_name,
                                dot_input=input_filename,
                                hyphy_batch_file=gy94bf1,
                                hyphy_result_file="gy94%s_%s.result" % (gene_id, ""))

        # gy model , part 1
        gy94p1_name = "%sgy94p1.bf" % gene_id
        gy94bf.set_partition(0, 60)
        gy94bf.write_batch_file(dot_aln=aln_name,
                                dot_input=input_filename,
                                hyphy_batch_file=gy94p1_name,
                                hyphy_result_file="gy94p1%s_%s.result" % (gene_id, ""))

        # gy model , part 2
        gy94p2_name = "%sgy94p2.bf" % gene_id
        gy94bf.set_partition(60, len_gene)
        gy94bf.write_batch_file(dot_aln=aln_name,
                                dot_input=input_filename,
                                hyphy_batch_file=gy94p2_name,
                                hyphy_result_file="gy94p2%s_%s.result" % (gene_id, ""))

        # nested model , full length
        nested_model.write_batch_file(dot_aln=aln_name,
                                      dot_input=input_filename,
                                      hyphy_result_file="nest_gy%s.result" % gene_id,
                                      hyphy_batch_file="%snest_gy.bf" % gene_id)

        # nested model , part 1
        nested_model.set_partition(0, 60)
        nested_model.write_batch_file(dot_aln=aln_name,
                                      dot_input=input_filename,
                                      hyphy_result_file="nest_gy%s_p1.result" % gene_id,
                                      hyphy_batch_file="%snest_gyp1.bf" % gene_id)

        # nested model , part 2
        nested_model.set_partition(60, len_gene)
        nested_model.write_batch_file(dot_aln=aln_name,
                                      dot_input=input_filename,
                                      hyphy_result_file="nest_gy%s_p2.result" % gene_id,
                                      hyphy_batch_file="%snest_gyp2.bf" % gene_id)

        # nt nested ,full
        bf_nt_nest.write_batch_file(dot_aln=aln_name,
                                    dot_input=input_filename,
                                    hyphy_result_file="%snt_nest.result" % gene_id,
                                    hyphy_batch_file="%snt_nest.bf" % gene_id)

        # nt nested , part1
        bf_nt_nest.set_partition(0, 60)
        bf_nt_nest.write_batch_file(dot_aln=aln_name,
                                    dot_input=input_filename,
                                    hyphy_result_file="%snt_nest_p1.result" % gene_id,
                                    hyphy_batch_file="%snt_nest_p1.bf" % gene_id)

        # nt nested , part2
        bf_nt_nest.set_partition(60, len_gene)
        bf_nt_nest.write_batch_file(dot_aln=aln_name,
                                    dot_input=input_filename,
                                    hyphy_result_file="%snt_nest_p2.result" % gene_id,
                                    hyphy_batch_file="%snt_nest_p2.bf" % gene_id)


        # gu model . full length
        gu_bf.write_batch_file(dot_aln=aln_name,
                               dot_input=input_filename,
                               hyphy_result_file="%sgu.result" % gene_id,
                               hyphy_batch_file="%sgu.bf" % gene_id)

        gu_bf.set_partition(0, 60)
        gu_bf.write_batch_file(dot_aln=aln_name,
                       dot_input=input_filename,
                       hyphy_result_file="%sgup1.result" % gene_id,
                       hyphy_batch_file="%sgup1.bf" % gene_id)

        gu_bf.set_partition(60 ,len_gene)
        gu_bf.write_batch_file(dot_aln=aln_name,
                       dot_input=input_filename,
                       hyphy_result_file="%sgup2.result" % gene_id,
                       hyphy_batch_file="%sgup2.bf" % gene_id)


def hyphy_traversal():
    # os.system("HYPHY %s >> %s" % (bf_filename, "%s_hyphy.log" % gene_id))
    bf_files = [file1 for file1 in os.listdir(os.curdir) if "bf" == file1.split(".")[-1]]

    for bf in bf_files:
        order = "HYPHY %s >> %s" % (bf, bf.split()[0] + "log")
        os.system(order)

if __name__ == "__main__":
    aln_folder_traversal("d:\\Work\\Custom\\compare0513")
    pass