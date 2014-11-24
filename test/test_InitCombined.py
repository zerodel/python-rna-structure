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



from slidingWindowDownStream import *
from InitCombined import *
import pyHYPHY.DataHYPHY as DH
import os
import os.path


def test_export_length():
    STEP_NUM = 21
    STEP_WIDTH_STEP = 15
    WINDOW_WIDTH_MIN = 30
    SITE_START_POINT = 0
    SITE_SHIFT_STEP = 15
    window_axis = [0,2]
    start_axis = range(STEP_NUM)

    result_path = "/Users/zerodel/WorkSpace/sliding300/gap"
    output_path = "/Users/zerodel/WorkSpace/sWAnalysis/300/gap/"
    with open(os.path.join(output_path, "seqLen.lst"), "w") as seq_len_table:
        
        seq_len_table.write("start_site/windowWidth\t"
                            + "\t".join(str(window_axis)) + "\n")
        for window_start_offset in range(STEP_NUM):
            len_seq_input = []
            for window_width_offset in window_axis:
                start_site = SITE_SHIFT_STEP*window_start_offset + SITE_START_POINT
                length_window = STEP_WIDTH_STEP*window_width_offset + WINDOW_WIDTH_MIN

                job_description = "w%ds%d" % (length_window, start_site)
                dot_input_file = os.path.join(result_path, "TIR%sN.input" % job_description)
                length_seq_in_input = input_file_info(dot_input_file)
                if not isinstance(length_seq_in_input, int):
                    print "not same length , exiting..."
                    return
                len_seq_input.append(str(length_seq_in_input))

            print len_seq_input
            seq_len_table.write("\t".join(len_seq_input) + "\n")


if __name__ == '__main__':
    test_export_length()
