# !/usr/bin/env python
# -*- coding : utf-8 -*-
# author :
# Readme:
# changelog :nothing
#

import re
import os
import logging
logging.basicConfig(filename=os.path.join(os.getcwd(), 'aln.log'), level=logging.DEBUG, filemode='w', format='%(asctime)s - %(levelname)s: %(message)s')


ecoli_pattern = "((SSO_XXXX:0.00250914321185456961,(SDY_XXXX:0.00512455264316328922,((SBO_XXXX:0.00255645533056974000,SXXXX:0.00297606412646914033):0.00048361360127168910,((GKPORF_BXXXX:0.05710356582382078439,(pluXXXX:0.22568960386950531749,SGXXXX:0.21100366205314127765):0.11856470733253710037):0.03145361970197192714,(STMXXXX:0.00180054347435949715,STYXXXX:0.00345947038924121399):0.03966386799217513220):0.03924531644843273076):0.00051331528259363401):0.00057744148977114274,bXXXX:0.00279499022522683872):0.0)"
ecoli_rooted_pattern = "(SSO_XXXX:0.00250914321185456961,(SDY_XXXX:0.00512455264316328922,((SBO_XXXX:0.00255645533056974000,SXXXX:0.00297606412646914033):0.00048361360127168910,((GKPORF_BXXXX:0.05710356582382078439,(pluXXXX:0.22568960386950531749,SGXXXX:0.21100366205314127765):0.11856470733253710037):0.03145361970197192714,(STMXXXX:0.00180054347435949715,STYXXXX:0.00345947038924121399):0.03966386799217513220):0.03924531644843273076):0.00051331528259363401):0.00057744148977114274,bXXXX:0.00279499022522683872):0.0"



def simple2(genes_share_alignment, pattern_2="(bXXXX,SSO_XXXX):0.0"):
    pattern_ecoli = "bXXXX"
    pattern_sso = "SSO_XXXX"

    pattern_all = [pattern_ecoli, pattern_sso]
    tree_content = pattern_2

    for index_p, item in enumerate(pattern_all):
        tree_content = re.sub(item, genes_share_alignment[index_p], tree_content)

    return tree_content

def ecoli(genes_share_alignment):
    """    make phylogenetic tree from zhoutong's ecoli data

    ((SSO_XXXX:0.00250914321185456961,(SDY_XXXX:0.00512455264316328922,((SBO_XXXX:0.00255645533056974000,SXXXX:0.00297606412646914033):0.00048361360127168910,((GKPORF_BXXXX:0.05710356582382078439,(pluXXXX:0.22568960386950531749,SGXXXX:0.21100366205314127765):0.11856470733253710037):0.03145361970197192714,(STMXXXX:0.00180054347435949715,STYXXXX:0.00345947038924121399):0.03966386799217513220):0.03924531644843273076):0.00051331528259363401):0.00057744148977114274,bXXXX:0.00279499022522683872):0.0)"
    """
    pattern_b = "bXXXX"         # pattern for ecoli genes
    pattern_sso = "SSO_XXXX"
    pattern_s = "SXXXX"
    pattern_sbo = "SBO_XXXX"
    pattern_sdy = "SDY_XXXX"
    pattern_gkporf = "GKPORF_BXXXX"
    pattern_stm = "STMXXXX"
    pattern_sty = "STYXXXX"
    pattern_plu = "pluXXXX"
    pattern_sg = "SGXXXX"

    pattern_all = [pattern_b, pattern_sso, pattern_s, pattern_sbo,
                   pattern_sdy, pattern_gkporf, pattern_stm, pattern_sty,
                   pattern_plu, pattern_sg]
    # ecoli_pattern is defined in the head of this file\
    # and now , use the rooted tree
    tree_content = ecoli_rooted_pattern
    if len(genes_share_alignment) != len(pattern_all):
        print "gene number not equal"
        logging.error("gene number not equal ----" + "number in aln:"
                      + str(len(genes_share_alignment)) + "provide in tree pattern:"
                      + str(len(pattern_all)))
        return
    else:

        for index_p, item in enumerate(pattern_all):
            # do something
            tree_content = re.sub(item, genes_share_alignment[index_p], tree_content)
    return tree_content

# todo : function ---> yeast : make tree for yeast
def yeast(genes_share_alignment):
    """

tree for yeast is ((Smik:0.06505810202920048191,(Skud:0.06824775484329334563,(Sbay:0.06449537146489706108,(Sklu:0.44406379138184037814,Scas:0.34852255969430995242):0.27361297838594561549):0.02621412500451008112):0.02467797569334097968):0.02448288696591136016,Spar:0.03007100476252094409,Scer:0.04122423044401372916):0.0;

    Args:

    Return:

    """
    ## just return the same string , it's ok

    return "((Smik:0.06505810202920048191,(Skud:0.06824775484329334563,(Sbay:0.06449537146489706108,(Sklu:0.44406379138184037814,Scas:0.34852255969430995242):0.27361297838594561549):0.02621412500451008112):0.02467797569334097968):0.02448288696591136016,Spar:0.03007100476252094409,Scer:0.04122423044401372916):0.0"



# todo : fly , make tree for fly
def fly(sequence_name, tree_pattern):
    pass


if __name__ == "__main__":
    pass