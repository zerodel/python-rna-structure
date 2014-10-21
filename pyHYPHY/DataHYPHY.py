"""
input file preparation of HYPHY

"""
__author__ = 'Zerodel_2'


class sequenceNotSameLength(Exception):
    """
    :raise when sequences have different length
    """
    pass

class notStr(Exception):
    pass

#  read dot-aln files and return the gene names in one aln file
def aln_reader(dot_aln_file):
    """
    read an aln file , and report genes in it .
    Args:
    1. dot_aln_file : file name of .aln file,which contains a alignment of several genes
    Return:
    1. genes: a list of strings , represent the genes
    """
    with open(dot_aln_file, "r") as aln:
        genes = [line.split()[0] for line in aln]

    return genes

def lines_to_input(aln_lines, input_file_path):
    """
    translate aln_lines for .aln file into a input file : input_file_path
    :param aln_lines:
    :param input_file_path:
    :return:
    """
    hyphy_input_content = [">%s\n%s\n" % (line.split()[0], line.split()[1])
                               for line in aln_lines]

    with open(input_file_path, "w") as Hyphy_input_writer:
        Hyphy_input_writer.writelines(hyphy_input_content)


def aln2input(dot_aln_file, hyphy_input_file):
    """translate dot-aln file to "input" file for HYPHY

    Args:
    dot_aln_file: file name of .aln alignment file
    hyphy_input_file : translated input file for hyphy process

    :Return:
    nothing
    """
    with open(dot_aln_file, 'r') as aln:

        hyphy_input_content = [">%s\n%s\n" % (line.split()[0], line.split()[1])
                               for line in aln]



    with open(hyphy_input_file, "w") as Hyphy_input_writer:
        Hyphy_input_writer.writelines(hyphy_input_content)


def aln2inputNogap(dot_aln_file, hyphy_input_file):
    with open(dot_aln_file, "r") as aln:
        contents = aln.readlines()
    gene_titles = [line.split()[0].strip() for line in contents]
    gene_sequences = [line.split()[-1].strip() for line in contents]

    check_sequence_matrix(gene_sequences)
    sequence_no_gap = remove_gaps_matrix(gene_sequences)

    with open(hyphy_input_file, "w") as write_hyphy:
        for index_i, gene_title in enumerate(gene_titles):
            write_hyphy.write(">" + gene_title.strip() + "\n")
            write_hyphy.write(sequence_no_gap[index_i].strip() + "\n")


def check_sequence_matrix(sequence_matrix):

    length_1 = len(sequence_matrix[0])

    for line in sequence_matrix:
        if not isinstance(line, str):
            raise notStr
        if not length_1 == len(line):
            raise sequenceNotSameLength


def remove_gaps_matrix(sequence_matrix):
    check_sequence_matrix(sequence_matrix)
    codon_sum = len(sequence_matrix[0])/3
    tmp_matrix_filter = [[] for line in sequence_matrix]

    # using a flipping window
    for codon_i in range(codon_sum):
        start_point = codon_i*3
        end_point = start_point + 3
        content_in_window = [line[start_point:end_point] for line in sequence_matrix]
        # if "-" exists , reject whole window
        if not "-" in "".join(content_in_window):
            # no "-"
            for index_line, line_filtered in enumerate(tmp_matrix_filter):
                line_filtered.append(content_in_window[index_line])

    sequence_matrix_filter = ["".join(line) for line in tmp_matrix_filter]
    return sequence_matrix_filter

def remove_gaps(sequence_string):
    """    remove gaps "---" in given sequence string ."""
    num_codon = len(sequence_string)/3
    seq_gap_free = []
    for codon_i in range(num_codon):
        codon_content = sequence_string[(codon_i*3):(codon_i*3 + 3)]
        if '-' in codon_content:
            continue
        seq_gap_free.append(codon_content)
    sequence_main = "".join(seq_gap_free)
    return sequence_main

def aln_info(aln_file):
    """read a .aln file and report gene names and length of sequence"""
    with open(aln_file, "r") as aln:
        infos = [(line.split()[0], len(line.split()[1])) for line in aln]

    genes = [info[0] for info in infos]
    lengths = [info[1] for info in infos]

    return genes, lengths


if __name__ == "__main__":
    pass