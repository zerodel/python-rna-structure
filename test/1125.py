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

# testing part in emacs
DEBUG_WITH_EMACS = True
if __name__ == "__main__":
    DEBUG_WITH_EMACS = False
if DEBUG_WITH_EMACS:
    print "ok"
    import os
    import os.path
    os.chdir("/Users/zerodel/WorkSpace/test")
    with open("0912psi.lst") as psi_reader:
        psi_content = psi_reader.readlines()

    with open("0912w.lst") as w_reader:
        w_content = w_reader.readlines()

        
    w_trans = ["%s\t%s\n" % (line[:5], line.split("=")[-1].strip()) for line in w_content]
    psi_trans = ["%s\t%s\n" % (line[:5], line.split("=")[-1].strip()) for line in psi_content]

    with open("psilst.txt", 'w') as psi_writer:
        psi_writer.write("geneid\tpsi\n")
        psi_writer.writelines(psi_trans)

    
    with open("wlst.txt", 'w') as psi_writer:
        psi_writer.write("geneid\tw\n")
        psi_writer.writelines(w_trans)

        
if __name__ == "__main__":
    pass
