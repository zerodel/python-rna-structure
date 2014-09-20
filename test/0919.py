# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#
import os
import os.path


def main1():
    path_2sp = "/home/zerodel/Workspace/merge0914/nest2psi.lst"
    path_10sp = "/home/zerodel/Workspace/merge0914/nest10psi.lst"
    path_factor_gene = "/home/zerodel/Workspace/merge0914/gene_factor_r.lst"

    with open(path_2sp, "r") as read2:
        content2 = read2.readlines()
        content2.pop(0)

    with open(path_10sp, "r") as read10:
        content10 = read10.readlines()
        content10.pop(0)

    genes2 = [line_single.split()[0] for line_single in content2]
    genes10 = [line_single.split()[0] for line_single in content10]

    with open(path_factor_gene, "w") as writer:
        writer.write("gene_id" + "\t" + "is_shared" + "\n")
        for gene_single in genes2:
            if gene_single in genes10:
                line_writen = gene_single + "\t" + "2" + "\n"
            else:
                line_writen = gene_single + "\t" + "1" + "\n"
            writer.write(line_writen)


if __name__ == "__main__":
    main1()