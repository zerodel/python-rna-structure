# -*- coding: utf-8 -*-
"""
Created on Sat Sep 20 20:49:38 2014

@author: zerodel
"""
import os
import os.path
import cPickle as pickle
import pyRNAsnp.pyRNAsnp as Rnasnp
import pyHYPHY.global_constants as gc

# Ecoli cpd was stored @ ~/Workspace/cpd_store_site
ecoli_cpd_folder = "/home/zerodel/Workspace/cpd_store_site"
codon_codon_list_file = "/home/zerodel/Workspace/codon_all.lst"
gene_id_file_path = "/home/zerodel/GitProjects/python-rna-structure/data/para/gu10w.txt"

def main1():
    # try to get a heat-map like graphic show .
    cpd_cpd_show_psi = Rnasnp.init_codon_pair_dict(gc.codon_list_64)
    with open(codon_codon_list_file, "r") as cclread:
        contents = cclread.read().split()

    psi_keys = [(ccl.split("_")[0], ccl.split("_")[-1]) for ccl in contents]

    for key in cpd_cpd_show_psi.keys():
        if key in psi_keys:
            cpd_cpd_show_psi[key].append(1)
        else:
            cpd_cpd_show_psi[key].append(0)
    with open("/home/zerodel/GitProjects/python-rna-structure/data/ecoli_psi.matrix", "w") as writer:

        codon_vector = gc.codon_list_64
        for codon1 in codon_vector:
            # build a line for lst file
            line_single = []
            for codon2 in codon_vector:
                line_single.append(str(cpd_cpd_show_psi[(codon1, codon2)][0]))

            writer.write("\t".join(line_single) + "\n")


def read_tab_get_gene_id(file_path):

    with open(file_path, "r") as fetch_txt:
        lines = fetch_txt.readlines()

    gene_id_order = [line[:5] for line in lines]
    return gene_id_order


def main2():
    # get a list of genes--> snp number
    #files_cpd = [cpd_file for cpd_file in os.listdir(ecoli_cpd_folder) if os.path.splitext(cpd_file)[1] == ".cpd"]
    #geneids = []
    geneids = read_tab_get_gene_id(gene_id_file_path)
    gene_lengths = []
    snp_nums = []
    for gene_id in geneids:
        #gene_id = cpd_file.split(".")[0]
        full_name = os.path.join(ecoli_cpd_folder, gene_id + ".cpd")
        fasta_name = os.path.join(ecoli_cpd_folder, gene_id + ".fasta")
            # load temp CPD
        with open(full_name, "r") as cpdf:
            cpd_temp = pickle.load(cpdf)
        num_snp = 0
        for key in cpd_temp.keys():
            num_snp += len(cpd_temp[key])
            # extend main CPD
        with open(fasta_name, "r") as fasta_reader:
            length_gene = len(fasta_reader.readlines()[-1])

        #geneids.append(gene_id)
        print full_name, "----- is done !"
        print fasta_name, "---- ready"
        snp_nums.append(num_snp)
        gene_lengths.append(length_gene)

    with open("/home/zerodel/GitProjects/python-rna-structure/data/ecoli_snp_num_gene_length.lst", "w") as writer:
        writer.write("geneid\tuseful_snp_num\tgene_length\n")
        for index_gi, gene_id in enumerate(geneids):
            line_2_write = "%s\t%s\t%s\n" % (gene_id, snp_nums[index_gi], gene_lengths[index_gi])
            writer.write(line_2_write)


if __name__ == "__main__":
    main2()