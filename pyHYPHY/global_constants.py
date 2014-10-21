codon_list = ['ATT',  'ATC',  'ATA',
              'CTT',  'CTC',  'CTA',  'CTG',  'TTA',  'TTG',
              'GTT',  'GTC',  'GTA',  'GTG',
              'TTT',  'TTC',
              'ATG',
              'TGT',  'TGC',
              'GCT',  'GCC',  'GCA',  'GCG',
              'GGT',  'GGC',  'GGA',  'GGG',
              'CCT',  'CCC',  'CCA',  'CCG',
              'ACT',  'ACC',  'ACA',  'ACG',
              'TCT',  'TCC',  'TCA',  'TCG',  'AGT',  'AGC',
              'TAT',  'TAC',
              'TGG',
              'CAA',  'CAG',
              'AAT',  'AAC',
              'CAT',  'CAC',
              'GAA',  'GAG',
              'GAT',  'GAC',
              'AAA',  'AAG',
              'CGT',  'CGC',  'CGA',  'CGG',  'AGA',  'AGG']

aa_codon_dict = {"I": ['ATT', 'ATC', 'ATA'],
                 "L": ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
                 "V": ['GTT', 'GTC', 'GTA', 'GTG'],
                 "F": ['TTT', 'TTC'],
                 "M": ['ATG'],
                 "C": ['TGT', 'TGC'],
                 "A": ['GCT', 'GCC', 'GCA', 'GCG'],
                 "G": ['GGT', 'GGC', 'GGA', 'GGG'],
                 "P": ['CCT', 'CCC', 'CCA', 'CCG'],
                 "T": ['ACT', 'ACC', 'ACA', 'ACG'],
                 "S": ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                 "Y": ['TAT', 'TAC'],
                 "W": ['TGG'],
                 "Q": ['CAA', 'CAG'],
                 "N": ['AAT', 'AAC'],
                 "H": ['CAT', 'CAC'],
                 "E": ['GAA', 'GAG'],
                 "D": ['GAT', 'GAC'],
                 "K": ['AAA', 'AAG'],
                 "R": ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']}

# synonymous codons
codon_aa_dict = {'ATT': 'I',  'ATC': 'I',  'ATA': 'I',
                 'CTT': 'L',  'CTC': 'L',  'CTA': 'L',  'CTG': 'L',  'TTA': 'L',  'TTG': 'L',
                 'GTT': 'V',  'GTC': 'V',  'GTA': 'V',  'GTG': 'V',
                 'TTT': 'F',  'TTC': 'F',
                 'ATG': 'M',
                 'TGT': 'C',  'TGC': 'C',
                 'GCT': 'A',  'GCC': 'A',  'GCA': 'A',  'GCG': 'A',
                 'GGT': 'G',  'GGC': 'G',  'GGA': 'G',  'GGG': 'G',
                 'CCT': 'P',  'CCC': 'P',  'CCA': 'P',  'CCG': 'P',
                 'ACT': 'T',  'ACC': 'T',  'ACA': 'T',  'ACG': 'T',
                 'TCT': 'S',  'TCC': 'S',  'TCA': 'S',  'TCG': 'S',  'AGT': 'S',  'AGC': 'S',
                 'TAT': 'Y',  'TAC': 'Y',
                 'TGG': 'W',
                 'CAA': 'Q',  'CAG': 'Q',
                 'AAT': 'N',  'AAC': 'N',
                 'CAT': 'H',  'CAC': 'H',
                 'GAA': 'E',  'GAG': 'E',
                 'GAT': 'D',  'GAC': 'D',
                 'AAA': 'K',  'AAG': 'K',
                 'CGT': 'R',  'CGC': 'R',  'CGA': 'R',  'CGG': 'R',  'AGA': 'R',  'AGG': 'R',
                 'TAA': '*',  'TAG': '*',  'TGA': '*',
                 '---': '-'}

codon_list_hypothesis = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA',
                         'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC',
                         'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG',
                         'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAC', 'TAT', 'TCA', 'TCC',
                         'TCG', 'TCT', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']


codon_list_64 = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA',
'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC',
'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG',
'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC',
'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']

codon_terminator = ('TAA', 'TAG', 'TGA')

dna_base_list = ('A', 'C', 'G', 'T')
rna_base_list = ('A', 'C', 'G', 'U')

is_purine = {"A":1,"G":1,"C":0,"T":0}