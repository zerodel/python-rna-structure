# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : Zerodel
# Readme:444
# functions for data analysis of RNAsnp outputs 

import sets
import os
import string
import global_constants as gc


try:
        import cPickle as pickle
except ImportError:
        import pickle


def rnasnp_extractor(rnasnp_file_name):
    """
    get the whole content of log file 
        argument:
        logfile: filename of .rnasnp file
        
        output 
        ----rnasnp_list, rnasnp_list_head----
        rnasnp_list: a  list contains the main content of log file
        rnasnp_list_head : a list contains the first row of log file ,  the head of table
    """
    with open(rnasnp_file_name, 'r') as file_rnasnp:
        rnasnp_list = file_rnasnp.readlines()
        rnasnp_list_head = rnasnp_list.pop(0)
    return rnasnp_list, rnasnp_list_head


def get_snp_list(rnasnp_list, snp_location=0):
    """
        get snps ( like  'A2959U' ) from .rnasnp file

        snps  =  get_snpList(loglist, snp_location = 0)
        argument ------------
        -loglist : output of logExtractor function 
        -snp_location :the location of snp data in a single line of log file
        output ------------
        -snps : a list of snp data ,  like "A423G"
    """
    # snps = []
    # for line in loglist:
    #     content_in_line = string.split(line)
    #     snptmp = content_in_line[snp_location]
    #     snps.append(snptmp)
    single_base_mutations = map(lambda rnasnp_line: rnasnp_line.split()[snp_location], rnasnp_list)
    return single_base_mutations


def codon_list_all(dna_base_list=gc.dna_base_list):
    """
    :param dna_base_list:
        Description:
        generate a codon list contain all possible codons ,  that is 4*4*4 = 64 , including 3 terminators
        Usage:
        codon_list =  codon_list_all(dna_base_list = dna_base_list):
        Arguments:
        dna_base_list:  a list contains the nucleinic acid bases,  ('A', 'C', 'G', 'T') by default
        Outputs:
        codon_list : a list contains all the possible triplets of amino acids ,  use the  bases in dna_base_list
        Changelog:
    """
    codon_list = []
    for c1 in dna_base_list:
        for c2 in dna_base_list:
            for c3 in dna_base_list:
                codon_triplet = ''.join([c1,  c2,  c3])
                codon_list.append(codon_triplet)
    return codon_list


def init_codon_pair_dict(list_triple_base=gc.codon_list):
    """

    :rtype : object
        Description:
        generate a codon -> codon  2D array , using dict
        codon -> codon relation , described like {("AAA", "AAC"):[]}
        Usage:
        CodonPairDict = CodonPairDict_init(list_triple_base = codon_list):
      
        Arguments:
        list_triple_base : list contains all preset codons,    output of codon_list_all by default
        OutPuts:
        a dict ,  elements like {("AAA", "AAC"):[]}
        Changelog:
    """
    ac_x = list_triple_base[:]
    ac_y = ac_x[:]
    dict_codon_pair = {}
    for xi in ac_x:
        for yi in ac_y:
            dict_codon_pair[(xi, yi)] = []

    return dict_codon_pair


def rnasnp_built_cpd(list_rnasnp, head_list_rnasnp, sequence_string, list_codon=gc.codon_list):
    """
    Description: 
    read a logfile of RNAsnp output ,  then built a CodonPairDict,
    
    Usage:
    CodonPairDict =  RNAsnp2CPD_M(list_rnasnp, head_list_rnasnp, sequence_string, list_codon = list_codon):
            
    arguments:
    list_rnasnp, head_list_rnasnp = data extracted from the Rnasnp output file ,  maybe XXX.rnasnp,
                        using logextractor method
    sequence_string = sequence string
    1. rebuild from log files
    2. read from fasta by default ,  using Bio.SeqIO.read,  one file ,  one sequence
             :"str(SeqIO.read(seqfile, 'fasta').seq)"
    Codonlist = list contains all the meaningful codon,  
       use output of codon_list_all by default
    OutPuts:
    CodonPairDict =  a dict ,          {("AAA", "AAC"):[xxx, xxxx, xxxx, xxxx]}
    
    Changelog:
    if the sequence_string contains "U",  that is ,  mRNA sequence,  error would happen
    add site information into these CPD

    2014-09-11: previous edition ignored the p-value ,
    use 0.1 as a threshold for p-value .

    """
    if "U" in sequence_string:
        trans_reverse = string.maketrans("U", "T")
        sequence_string = sequence_string.translate(trans_reverse)

    dict_codon_pair = init_codon_pair_dict(list_codon)
    data_position = get_position_data_in_record(head_list_rnasnp, "d_max")
    # define a p-value position and a threshold for p-value
    p_value_threshold = 0.1
    p_value_pos = -1

    for record in list_rnasnp:
        list_row = record.split()
        snp_record = list_row[0]
        snp_site = int(snp_record[1:-1])
        data_record = list_row[data_position]
        p_value = float(list_row[p_value_pos])
        codon_origin = get_snp_codon_origin(snp_record, sequence_string)
        codon_target = get_snp_codon_target(snp_record, sequence_string)
        if "N" in codon_origin or "N" in codon_target:
                continue
        if codon_origin in gc.codon_terminator or codon_target in gc.codon_terminator:
                continue
        elif len(data_record) > 0 and p_value < p_value_threshold:
            dict_codon_pair[(codon_origin, codon_target)].append((data_record, snp_site))  # worth mention here
    return dict_codon_pair


def get_position_data_in_record(head_list_rnasnp, data_name='d_max'):
    """
    return required data position (column ) in a single row of log files (.rnasnp files)
    head_list_rnasnp is default like this
    "SNP	w	Slen	GC	interval	d	pvalue1	ewin	interval	d_max	pvalue2"
    """
    list_head_list_rnasnp = head_list_rnasnp.split("\t")
    return list_head_list_rnasnp.index(data_name)


def snp_codon_position_translate(snp_record, string_sequence):
    """
        Description: 
        get the codon beginning site in the string_sequence
        
        Usage:
        Position_begin_codon = snp_codon_position_translate(snp_record, string_sequence)
        
        arguments:
        snp_record : like "A3234G" in rnasnp output files,  
        string_sequence: sequence string read from .fasta files
        OutPuts:
        Position_begin_codon : an integer marking the beginning of the codon where the snp occurred
        
        Changelog:
        
    """
    position_snp = int(snp_record[1:-1]) - 1
    index_codon = position_snp/3
    position_begin_codon = index_codon * 3
    position_end_codon = position_begin_codon + 3
    if position_end_codon > len(string_sequence):
        position_begin_codon = len(string_sequence) - 3
    return position_begin_codon


def get_snp_codon_origin(snp_record,  string_sequence):
    """
        Description: 
        get the origin  codon of the very snp_record in the string_sequence
        attention: the position in snp_record begins at 1 ,  not 0 (as in Python)

        
        Usage:
        Codon_Origin = get_snp_codon_origin(snp_record, string_sequence):
        
        arguments:
        snp_record : like "A3234G" in rnasnp output files,  
        string_sequence: sequence string read from fasta files
        
        OutPuts:
        Codon_Origin : the triplet of the origin codon in snp mutation
        
        Changelog:
        
    """
    position_begin_codon = snp_codon_position_translate(snp_record, string_sequence)
    codon_origin = string_sequence[position_begin_codon: position_begin_codon + 3]
    return codon_origin


def get_snp_codon_target(snp_record, string_sequence):
    """
        Description: 
        get the target codon of the very snp_record in the string_sequence
        attention: the position in snp_record begins at 1 ,  not 0 (as in Python)

        
        Usage:
        Codon_Target = get_snp_codon_origin(snp_record, string_sequence):
        
        arguments:
        snp_record : like "A3234G" in rnasnp output files,  
        string_sequence: sequence string read from fasta files
        
        OutPuts:
        Codon_Target : the triplet of the 'target' codon in snp mutation
        
        Changelog:
        
    """
    position_snp = int(snp_record[1:-1]) - 1
    snp_in_codon = position_snp % 3
    codon_in_seq = get_snp_codon_origin(snp_record, string_sequence)
    codon_listed = list(codon_in_seq)
    nt_snp = snp_record[-1]
    if 'U' == nt_snp:
        nt_snp = 'T'
    codon_listed[snp_in_codon] = nt_snp
    codon_target = ''.join(codon_listed)
    return codon_target


def print_codon(dict_data_codon_pair, codon_str, filename_data, site_num=-1):
    """
        Description: 
           write all the data(in CPD) of the specified codon (with the name of Codon) to file( filename_data)
        
        Usage:
        filename_data = Print_Codon(CPD, codon_str, filename_data)
        
        arguments:
        -CPD =  data dict ,  output of RNAsnp2CPD
        -codon_str : the chosen codon
        -filename_data : filename of output data tab file
        Outputs:
        -filename_data = filename contains the tab data 
        -site_num: decide how many snps should be processed downstream , if set to -1 , all snps will be used 
        Changelog:
        
    """
    list_codon_string = list(codon_str)
    # check codon substitution direction .
    if list_codon_string[0] == '_' and list_codon_string[-1] != '_':
        position_in_tuple = 1
        string_codon_triplet = ''.join(list_codon_string[1:])
    elif list_codon_string[-1] == '_' and list_codon_string[0] != '_':
        position_in_tuple = 0
        string_codon_triplet = ''.join(list_codon_string[:-1])
    else:
        print "error: wrong input codon string .. [ _AAA or AAA_ is right]"
        return None
    print string_codon_triplet, '\n', '========================='

    data_write_to_file = []

    for Codon_key in dict_data_codon_pair.keys():
        if Codon_key[position_in_tuple] == string_codon_triplet and len(dict_data_codon_pair[Codon_key]) > 0:
            if -1 == site_num:
                data_content = [x[0] for x in dict_data_codon_pair[Codon_key]]
            elif site_num > 0:
                data_content = [x[0] for x in dict_data_codon_pair[Codon_key] if x[1] < site_num]
                codon_content = str(Codon_key[0]) + "_" + str(Codon_key[1])
                file_line = codon_content + " " + " ".join(data_content)
                data_write_to_file.append(file_line)

    file_content = '\n'.join(data_write_to_file)
    try:
        file_handler = open(filename_data, 'w+')
    except Exception, e:
        print "Error in opening file "
        return str(e)
    else:
        file_handler.writelines(file_content)
    finally:
        file_handler.close()
    return filename_data


def rebuild_seq(rnasnp_file):
    """
        Description:
         - rebuild a sequence using a RNAsnp output file .
        Arguments:
        - `rnasnp_file`: filename or path to the RNAsnp file
        
        Output:
        -'seq_String' : the rebuilt sequence 
    """
        # read snps 
    rnasnp_list, rnasnp_list_head = rnasnp_extractor(rnasnp_file)
    snps = get_snp_list(rnasnp_list)

#    seq_nt = []
#    seq_position = []
#    for snp in snps:
#        seq_nt.append(snp[0])
#        seq_position.append(int(snp[1:-1]))
    seq_nt = [snp[0] for snp in snps]
    seq_position = [int(snp[1:-1]) for snp in snps]
    seq_vector = zip(seq_position, seq_nt)
    seq_unique = list(sets.Set(seq_vector))
        # bubble sort
        # b=zip(a.keys(),  a.values())    # generate a list of tuples
        # sorted(b,  key=lambda item:item[1])  
        # rebuild the seq
        # seqR = [pointer[1] for pointer in seq_unique]
    sort_tuple_list(seq_unique)  # in site bubble sort for  list made of  tuples
    seqR = []
        
    for index_nt in range(len(seq_unique)):
        if seq_unique[index_nt][0] == index_nt + 1:
            seqR.append(seq_unique[index_nt][1])
        else:
            seqR.append("N")
    seq_string = "".join(seqR)
    return seq_string


def sort_tuple_list(seq_unique, pt=0):
    """
    Description:
    in-site bubble sort for list made of tuples
    Usage:
    suppose original list like:[(1, A), (3, B), (2, D)....] ,  sorted by number in (1, A) ,  pt should be 0
    
    Arguments:
    - `seq_unique`: list of tuple(list),  elements in list should be in the same format
    - 'pt': the index position in tuple 
    Outputs:
    - None 
    """
    for ii in range(len(seq_unique)):
        for jj in range(len(seq_unique)):
            if ii < jj and seq_unique[ii][pt] > seq_unique[jj][pt]:
                seq_unique[ii], seq_unique[jj] = seq_unique[jj], seq_unique[ii]
    
    return None


def snp_dir_list(rnasnp_output_path):
    """
    accept a file or a folder name , and return a list of .rnasnp files and folder name of .rnasnp files
    Arguments:
    - `rnasnp_output_path`:
    """
    if os.path.isfile(rnasnp_output_path):  # check whether the parameter given is a file
        rnasnp_folder_name, rnasnp_filename = os.path.split(rnasnp_output_path)
        if not "" == rnasnp_folder_name:  # check whether the file name is only a base name
            list_all_files = os.listdir(rnasnp_folder_name)
        else:
            list_all_files = os.listdir(".")
    elif os.path.isdir(rnasnp_output_path):  # if given a folder
        list_all_files = os.listdir(rnasnp_output_path)
        rnasnp_folder_name = rnasnp_output_path
    else:
        print "please input correct name for  file or folder"
        return None
    rnasnp_files = [filename for filename in list_all_files if os.path.splitext(filename)[1] == ".rnasnp"]
    return rnasnp_files, rnasnp_folder_name


def snp_dir_traversal(rnasnp_files, rnasnp_folder_name, output_folder="", is_display_progress=False, is_force_update=True):
    """
    do a traversal over .rnasnp files' folder. to each .rnasnp file, build a cpd structure and dump it to .cpd file.

    Arguments:
    - `rnasnp_folder_name`:
    - `rnasnp_files`:
    """
    num_seqn = 0
    num_cpd_done = 0
    num_cpd_all = len(rnasnp_files)
    for single_rnasnp_file in rnasnp_files:
            # single file operation
        whole_rnasnp_file_name = os.path.join(rnasnp_folder_name, single_rnasnp_file)
        rnasnp_basename_no_suffix = os.path.splitext(single_rnasnp_file)[0]
        seq_file_name = os.path.join(output_folder, "".join([rnasnp_basename_no_suffix, ".fasta"]))
        cpd_file_name = os.path.join(output_folder, "".join([rnasnp_basename_no_suffix, ".cpd"]))

        rnasnp_list, rnasnp_list_head = rnasnp_extractor(whole_rnasnp_file_name)

        if not os.path.exists(seq_file_name) or is_force_update:  # if .fasta file not exist or forced to update

            seq_string = rebuild_seq(whole_rnasnp_file_name)
            
            if "N" in seq_string:
                print "===>>", single_rnasnp_file, string.index(seq_string, "N")
                num_seqn += 1
         
            with open(seq_file_name, "w") as seq_file:
                if seq_file:
                    fasta_head = "> rebuild sequence - "
                    seq_file.write(fasta_head)
                    seq_file.write(seq_string)
            
        if not os.path.exists(cpd_file_name) or is_force_update:
            with open(seq_file_name, "r") as seq_reader:
                if seq_reader:
                    seq_string = seq_reader.readlines()[1]
                else:
                    print "error: sequence file error"
                    return None
            codon_pair_dict_single_file = rnasnp_built_cpd(rnasnp_list, rnasnp_list_head, seq_string)
                    
            with open(cpd_file_name, "w") as CPD_file:
                if CPD_file:
                    pickle.dump(codon_pair_dict_single_file, CPD_file)
                    num_cpd_done += 1
                    if is_display_progress:
                        print str(num_cpd_done), "files Done ,  Total number is ", str(num_cpd_all)
                        print single_rnasnp_file, "=====",  seq_file_name, "----", cpd_file_name
    print "total sequences number that contains N is ", num_seqn
    return None

def cpd_dir_traversal(codon_str, cpd_dir, output_dir, site_num=-1):
    """
    collect "codon_name" data from "cpd_dir" and output to "output_dir"
    Arguments:
    - `codon_str`:
    - `cpd_dir`:
    - `output_dir`:
    """
    # get whole list of cpd_dir
    files_in_cpd_dir = os.listdir(cpd_dir)
    # get the list of .cpd files
    files_cpd = [cpd_file for cpd_file in files_in_cpd_dir if os.path.splitext(cpd_file)[1] == ".cpd"]
    # set the codon name,  get the Codon_str and position_in_tuple
    list_codon_triplet = list(codon_str)
    if list_codon_triplet[0] == '_' and list_codon_triplet[-1] != '_':  # check codon substitution direction
        position_in_tuple = 1
        codon_str = ''.join(list_codon_triplet[1:])
    elif list_codon_triplet[-1] == '_' and list_codon_triplet[0] != '_':
        position_in_tuple = 0
        codon_str = ''.join(list_codon_triplet[:-1])
    else:
        print "error: wrong codon string"
        return None
    # initialize the main CPD dict
    cpd_main = init_codon_pair_dict(list_triple_base=gc.codon_list)
    # process each .cpd file ,  load temp CPD dict,  and extend main CPD with temp CPD item 
    for cpd_file in files_cpd:
        full_name = os.path.join(cpd_dir, cpd_file)
            # load temp CPD
        with open(full_name, "r") as cpdf:
            cpd_temp = pickle.load(cpdf)
            # extend main CPD
        for Codon_key in cpd_temp.keys():
            if Codon_key[position_in_tuple] == codon_str and len(cpd_temp[Codon_key]) > 0:
                cpd_main[Codon_key].extend(cpd_temp[Codon_key])
    # write main CPD to the output file .
    filename_data = os.path.join(output_dir, codon_str + ".lst")
    filename_data = print_codon(cpd_main, codon_str, filename_data, site_num)
    return cpd_main


def get_whole_cpd(cpd_dir):
    """
    collect all data from all .cpd files in "cpd_dir" and put into one cpd dict structure
    
    Arguments:
    - `cpd_dir`:
    """
    # get whole list of cpd_dir
    list_files_in_cpd_dir = os.listdir(cpd_dir)
    # get the list of .cpd files
    list_cpd_files = [cpd_file for cpd_file in list_files_in_cpd_dir if os.path.splitext(cpd_file)[1] == ".cpd"]
    cpd_main = init_codon_pair_dict(list_triple_base=gc.codon_list)
    # process each .cpd file ,  load temp CPD dict,  and extend main CPD with temp CPD item 
    for cpd_file in list_cpd_files:
        full_name = os.path.join(cpd_dir, cpd_file)
            # load temp CPD
        with open(full_name, "r") as cpdf:
                cpd_temp = pickle.load(cpdf)
            # extend main CPD
        for Codon_key in cpd_temp.keys():
            if len(cpd_temp[Codon_key]) > 0:
                cpd_main[Codon_key].extend(cpd_temp[Codon_key])
    return cpd_main


def cpd_codon_quantile_stat(cpd_quantile, list_codon=gc.codon_list, is_origin=False):
    """
    """
    quantile_mean = {}
    if is_origin:
        position_in_key_tuple = 0
    else:
        position_in_key_tuple = 1

    for codon in list_codon:
        sum_quantile = 0
        num_codon_pair = 0
        for cpd_key in cpd_quantile.keys():
            if isinstance(cpd_quantile[cpd_key], float) and codon == cpd_key[position_in_key_tuple]:
                sum_quantile += cpd_quantile[cpd_key]
                num_codon_pair += 1
        quantile_mean.setdefault(codon, sum_quantile/float(num_codon_pair))
    return quantile_mean


def build_state_cpd(dict_cpd, list_codon=gc.codon_list, larger_quantile_cc=[]):
    """
    TODO: need more modification , not that easy
     Arguments:
     - `CPD`:
    """
    cpd_state = init_codon_pair_dict(list_codon)
    if 0 == len(larger_quantile_cc):
        cpd_quantile = cpd_codon_quantile(dict_cpd)
        quantile_mean = cpd_codon_quantile_stat(cpd_quantile, gc.codon_list, False)

    for key in cpd_state.keys():
        codon_origin = key[0]
        codon_target = key[1]
        if key[0] == key[1]:
            cpd_state[key] = "*"
        elif len(dict_cpd[key]) == 0:
            cpd_state[key] = "0"
        elif gc.codon_aa_dict[codon_origin] != gc.codon_aa_dict[codon_target]:
            nt_pair = compare_codon_pair(key)
            nt_pair.sort()
            cpd_state[key] = "".join(["a", "".join(nt_pair), "*", "w"])
        elif gc.codon_aa_dict[codon_origin] == gc.codon_aa_dict[codon_target]:
            nt_diff = compare_codon_pair(key)
            nt_diff.sort()
            nt_pair = "".join(nt_diff)
            if 0 == len(larger_quantile_cc):
                if cpd_quantile[key] > quantile_mean[codon_origin]:
                    cpd_state[key] = "".join(["a", nt_pair, "*", "psi"])
                else:
                    cpd_state[key] = "".join(["a", nt_pair])
            else:
                if "_".join([codon_origin, codon_target]) in larger_quantile_cc:
                    cpd_state[key] = "".join(["a", nt_pair, "*", "psi"])
                    # describe the psi parameter 
                    print cpd_state[key], "==", nt_diff
                else:
                    cpd_state[key] = "".join(["a", nt_pair])

    return cpd_state


def diff_codon(codon1, codon2):
    """
    report how many differences between two codons
    """
    length_codon = 3
    diff_occur = []
    for index_in_codon in range(length_codon):
        if codon1[index_in_codon] == codon2[index_in_codon]:
            continue 
        else:
            diff_occur.append(index_in_codon)
    
    return diff_occur


def compare_codon_pair(codon_pair):
    """
    CodonPair: a tuple contain 2 codon ,and one difference
    example:("AAA","AAC") or ("AAC","AAA")
    all return ["A","C"]
    """
    codon1 = codon_pair[0]
    codon2 = codon_pair[1]

    if codon1 == codon2:
        return []

    nt_diff = []
    diff_occur = diff_codon(codon1, codon2)

    if 1 == len(diff_occur):
        diff = diff_occur[0]
        nt_diff.append(codon1[diff])
        nt_diff.append(codon2[diff])
        return nt_diff
    else:
        return None


def cpd_codon_quantile(dict_cpd):
    """
    
    Arguments:
    - `CPD`:
    - `codon_str`:
    - `filename_data`:
    """
    cpd_quantile = init_codon_pair_dict()
    for Codon_key in dict_cpd.keys():
        if len(dict_cpd[Codon_key]) > 0:
            data_content_float = [float(list_key_content) for list_key_content in dict_cpd[Codon_key]]
            quantile_float = quantile_of_list(data_content_float, 95)
            cpd_quantile[Codon_key] = quantile_float
    return cpd_quantile


def quantile_of_list(float_list, percents=95):
    """
    Arguments:
    - `float_list`:
    - 'percents': an integer
    """
    list_sorted = sorted(float_list)
    pointerx = percents*len(list_sorted)/100
    return list_sorted[pointerx]


def get_large_codon_group(filename):
    """ get codon - codon pairs for a file containing "AAC_TTT" like strings
    """
    larger_c2c_tuples = []

    with open(filename, "r+") as fileop:
        for line in fileop:
            larger_c2c_tuples.extend(line.split())
    return larger_c2c_tuples


def is_synonymous(codon_codon_pair):
    codon1 = codon_codon_pair[0:3]
    codon2 = codon_codon_pair[4:]
    return gc.codon_aa_dict[codon1] == gc.codon_aa_dict[codon2]


def mutation_one_site(single_base, position):
    """
    produce list of single site mutations with length 3.
    :param single_base: single nucleic acid base, within {A,G,T,C}
    :param position: the position of mutation in the sequence
    :return: a list of mutation in one site, like A3G , 3 indicate the position
    """
    dna_nt_4 = ('A', 'C', 'G', 'T')
    # output_list = []
    # for nt in dna_nt_4:
    #     if not single_base == nt:
    #         output_list.append("".join([single_base, str(position), nt]))

    # a "smarter" way to solve it
    output_list = ["".join([single_base, str(position), nt]) for nt in dna_nt_4 if not single_base == nt]
    return output_list



#####---------------------------------------------------------------------------
####################################################################
if __name__ == "__main__":
    pass
    # parent_path = "/home/zerodel/workspace/ZRna/"
    # source_path = "/home/zerodel/workspace/localdata/ecoli"
    # output_path = parent_path + "cpd_store_site"
    # single_codon_file_path = parent_path + "50single_codon_lst"
    # # rnasnp_files, snppathname = snp_dir_list(source_path)
    # #step 1 , set up  .cpd files
    # # snp_dir_traversal(rnasnp_files, snppathname, output_path, single_codon_file_path)
    # #step 2 , print each codon lst file
    # # for codon in codon_list_hypothesis:
    # #     codon_with_direction = codon + "_"
    # #     print "--->", codon_with_direction
    # #     cpd_dir_traversal(codon_with_direction, output_path, single_codon_file_path, 60)
    #
    # # # R kmeans part ******
    # # print "kmeans of cluster"
    # # kmeans_script_path = os.path.join(parent_path + "src","kmeans.R")
    # cc_significant_file_path = os.path.join(parent_path + "src", "codon20.lst")
    # # cmd_line = "Rscript " + kmeans_script_path + " " +  single_codon_file_path + " > " + cc_significant_file_path
    # # os.system(cmd_line)
    #
    # # read the result of R process
    # cc_significant_raw = get_large_codon_group(cc_significant_file_path)
    # cc_significant = [x for x in cc_significant_raw if is_synonymous(x)]
    # print len(cc_significant)
    #
    # CPD_site = get_whole_cpd(output_path)
    # CPDState = build_state_cpd(CPD_site, gc.codon_list_hypothesis, cc_significant)
    # # print_state_cpd(CPDState,"60.def",codon_list_hypothesis)
    #
    # # if not os.path.exists("CPD_site.db"):
    # #     with open("CPD_site.db", "w") as writecpd:
    # #         pickle.dump(CPD_site, writecpd)
    # # else:
    # #     with open("CPD_site.db", "r") as fcpd:
    # #         CPD_site = pickle.load(fcpd)
    # # print "----"
    # # for key in CPDState.keys():
    # #     if "psi" in CPDState[key]:
    # #         if compareCodonPair(key)[1] == "T":
    # #             print key