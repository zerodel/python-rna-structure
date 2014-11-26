# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:
#
__author__ = 'Zerodel'



def AddSysPath(new_path):
    """ AddSysPath(new_path): adds a "directory" to Python's sys.path
    Does not add the directory if it does not exist or if it's already on
    sys.path. Returns 1 if OK, -1 if new_path does not exist, 0 if it was
    already on sys.path.
    """
    import sys, os
    # Avoid adding nonexistent paths
    if not os.path.exists(new_path): return -1
    # Standardize the path.  Windows is case-insensitive, so lowercase
    # for definiteness if we are on Windows.
    new_path = os.path.abspath(new_path)
    if sys.platform == 'win32':
        new_path = new_path.lower( )
    # Check against all currently available paths
    for x in sys.path:
        x = os.path.abspath(x)
        if sys.platform == 'win32':
            x = x.lower( )
        if new_path in (x, x + os.sep):
            return 0
    sys.path.append(new_path)
    # if you want the new_path to take precedence over existing
    # directories already in sys.path, instead of appending, use:
    # sys.path.insert(0, new_path)
    return 1

AddSysPath("..")


import pyRNAsnp.pyRNAsnp as pyr
import pyHYPHY.global_constants as gc

# testing part in emacs
DEBUG_WITH_EMACS = True
if __name__ == "__main__":
    DEBUG_WITH_EMACS = False
if DEBUG_WITH_EMACS:
    print "hello world"
    cpd_path = "/Users/zerodel/WorkSpace/data/cpd_full"
    whole_cpd = pyr.get_whole_cpd(cpd_path)
    length_cc = [(key, len(whole_cpd[key])) for key in whole_cpd.keys() if len(whole_cpd[key])]
    length_cc = sorted(length_cc)
    codon_order = gc.codon_list_hypothesis
    with open("length_cpd.lst", "w") as writek:
        writek.write("\t".join(codon_order) +"\n")
        for c1 in codon_order:
            for c2 in codon_order:
                writek.write(str(len(whole_cpd[(c1,c2)])))
                writek.write("\t")
            wri
if __name__ == "__main__":
    pass
