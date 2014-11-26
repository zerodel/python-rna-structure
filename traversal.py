import os

def hyphy_traversal():
    # os.system("HYPHY %s >> %s" % (bf_filename, "%s_hyphy.log" % gene_id))
    bf_files = [file1
             for file1 in os.listdir(os.curdir)
             if "bf" == file1.split(".")[-1]]
    print len(bf_files)

    for bf in bf_files:
        order = "HYPHY %s | tee %s" % (bf , bf.split()[0] + "log")
        print order
        os.system(order)


if __name__ == "__main__":
    hyphy_traversal()
