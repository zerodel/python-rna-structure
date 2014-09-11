# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author :
# Readme:
# changelog :nothing
#
__author__ = 'Zerodel_2'

import re
import global_constants as gc

import logging
logging.basicConfig(filename="mdlHYPHY.log", filemode="a", level=logging.DEBUG)

class ModelHYPHY(object):
    """produce a .mdl file for HYPHY.
    have some switches for
    1. use rna structure information or not
    2. use detailed nt names or just mark syn or nonsyn
    """

    def __init__(self, lst_2_tuple_file_name="", is_rna_structure=False, is_nt_specific=False, special_tuple_list=None):

        """        constructor for mdlHYPHY        """
        self.model_content = []
        self.larger_codon_tuples = []
        self.is_nt_detailed = is_nt_specific
        self.is_rna_structure = is_rna_structure
        if not None == special_tuple_list:
            self.larger_codon_tuples = special_tuple_list
        elif not "" == lst_2_tuple_file_name:
            self.larger_codon_tuples = self.read_lst(lst_2_tuple_file_name)

        self._custom_model_definition()

    @staticmethod
    def read_lst(lst_2_tuple_file_name):
        total_codon_tuples_string = []
        with open(lst_2_tuple_file_name, "r") as lst_reader:
            for text_line in lst_reader:
                total_codon_tuples_string.extend(text_line.split())

        total_codon_tuples = [(tuple_codon[:3], tuple_codon[4:]) for tuple_codon in total_codon_tuples_string]

        logging.debug("total number of codon tuples is %d" % len(total_codon_tuples))
        return total_codon_tuples

    @staticmethod
    def _diff_codon_tuple(codon1, codon2):
        """        compare two codons.
        if the same , return None
        else return difference        """
        tt_codon = zip(codon1, codon2)
        difference_codon = [tt for tt in tt_codon if tt[0] != tt[1]]
        return difference_codon

    def check_codons(self, codon_origin="", codon_target="", codon_tuple=None):
        """
        check codon tuples codon_origin to codon_target.
        if origin and target are identical . return a bool value True
        if origin and target have more than 1 difference in their 3 nts .
             showing how many nt they differs

        if origin and target share two identical nts , this function will return a tuple
        (is_transversion, is_sysnonmous, is_larger, nt_diff)
        is_transversion , shows whether the biologic mutation is a transversion .
        is_sysnonmous , shows whether the mutation is sysnonmous
        is_larger , shows whether the mutation result in change of  rna structure
        nt_diff , returns the different nts .
        """
        if not None == codon_tuple:
            codon_origin = codon_tuple[0]
            codon_target = codon_tuple[1]

        if codon_origin == codon_target:
            identical_codons = True
            return identical_codons

        diff_in_tuple = self._diff_codon_tuple(codon_origin, codon_target)
        logging.debug("diff is %s" % str(diff_in_tuple))
        if len(diff_in_tuple) > 1:  # two codon have more than one difference
            return len(diff_in_tuple)
        else:
                diff_nt = list(diff_in_tuple[0])
                diff_nt.sort()

                is_transversion = gc.is_purine[diff_nt[0]] ^ gc.is_purine[diff_nt[1]]
                is_sysnonmous = gc.codon_aa_dict[codon_origin] == gc.codon_aa_dict[codon_target]
                is_larger = (codon_origin, codon_target) in self.larger_codon_tuples
                nt_diff = "".join(diff_nt)
                return is_transversion, is_sysnonmous, is_larger, nt_diff

    def _custom_model_definition(self):
        """build a HYPHY model definition following  syntax.
                define a element of matrix each line
        """
        codon_axis = gc.codon_list_hypothesis
        logging.debug("length of codon_list is %d" % len(codon_axis))

        definition_lines = []
        shift_x = 0
        shift_y = 0
        for index_origin in range(len(codon_axis)):
            #index_origin -= 1
            codon_origin = codon_axis[index_origin]
            if codon_origin in gc.codon_terminator:
                shift_x += 1
                continue

            for index_target in range(len(codon_axis)):
            #    index_target -= 1
                codon_target = codon_axis[index_target]

                if codon_target in gc.codon_terminator:
                    shift_y += 1
                    continue

                # here comes (codon1, codon2) processing
                codon_tuple_result = self.check_codons(codon_origin, codon_target)
                if bool == type(codon_tuple_result):  # two codons are identical
                    continue

                elif int == type(codon_tuple_result):  # two codons differs more than 2 nts
                    continue

                else:
                    is_transversion, is_sysnonmous, is_larger, nt_diff = codon_tuple_result
                    substitution_description = []
                    if self.is_nt_detailed:
                        substitution_description.append(nt_diff)
                    else:
                        if is_sysnonmous:
                            substitution_description.append("syn")
                        else:
                            substitution_description.append("nonsyn")

                    if not is_transversion:
                        substitution_description.append("k")

                    if self.is_rna_structure and is_larger:
                        substitution_description.append("psi")

                    model_description_line = "\tModelMatrixName[%d][%d] := %s ;" % \
                                             (index_origin - shift_x, index_target - shift_y, "*".join(substitution_description))

                    definition_lines.append(model_description_line)

        self.model_content = definition_lines
        return definition_lines

    def export_model(self, template_file_name="", target_mdl_name="", matrix_content=None):

        with open(template_file_name, "r") as template_reader:
            mdl_template = template_reader.read()

        if None == matrix_content and [] != self.model_content:  # use class variable
            matrix_content = self.model_content

        to_insert = "\n".join(matrix_content)
        mdl_template = re.sub("CUSTOM_MODEL_DEFINITION", to_insert, mdl_template)

        with open(target_mdl_name, "w") as mdl_writer:
            mdl_writer.write(mdl_template)

        return to_insert

if __name__ == "__main__":
    pass
