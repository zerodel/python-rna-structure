"""
    classes for making HYPHY .bf file

usage:
        gy94mdl = "/usr/lib/hyphy/TemplateBatchFiles/TemplateModels/GY94.ModelH"
        own_nest_mdl = "gy_nest_60.ModelH"
        bf_maker = bfH.HYPHYBatchFile(species_name="ecoli",
                                      model_file=own_nest_mdl,
                                      bf_template_file="../pyHYPHY/batch_template.bf")

        bf_maker.write_batch_file(dot_input="b0099.input",
                                  dot_aln="b0099.aln",
                                  hyphy_result_file="nested_gy.result",
                                  hyphy_batch_file="nested_gy94.bf")

"""
__author__ = 'Zerodel_2'

import re
import os.path
import pyHYPHY.SpeciesSpecificTree as SSTree
import pyHYPHY.DataHYPHY as pHdata




class HYPHYBatchFile(object):
    """
    main interface of bf handler
    """

    def __init__(self, species_name="", model_file="", bf_template_file=""):
        """
        constructor for bfHandler
        """
        self.species = species_name

        self.mdl_file = model_file
        self.bf_template = bf_template_file

        self.f_input = "INPUT_FILE"
        self.f_output = "OUTPUT_RESULT_FILE"
        self.f_mdl = "CUSTOM_MODEL_FILE"
        self.f_tree = "TREE_NEWICK_STRING"
        self.f_partition = "PARTITION_DEFINITION"
        self.f_matrix_name = "MATRIX_NAME_FOR_MODEL"
        self.matrix_name = self.get_matrix_name(self.mdl_file)

        self.use_given_tree = False
        self.tree_definition_external = ""

        self.partition = (0, 0)

        if "" != self.bf_template:
            with open(name=self.bf_template, mode="r") as bf_template_reader:
                self.batch_content = bf_template_reader.read()
        else:
            self.batch_content = ""

    def set_tree_from_outside(self, newick_string):
        self.use_given_tree = True
        self.tree_definition_external = newick_string

    # get a matrix name from your mdl definition
    def get_matrix_name(self, mdl_file):
        """ read a .mdl file and return ONE matrix name or a matrix name list
        now , just one
        :param mdl_file:
        :return:
        """
        with open(mdl_file, "r") as mdl_reader:
            matrix_name = [line.split("=")[0].strip() for line in mdl_reader if "=" in line]

        return matrix_name

    def set_partition(self, begin=0, end=0):
        """ set the partition range of sequence """
        self.partition = (begin, end)

    def build_tree(self, genes_share_one_alignment):
        """make phylogenetic tree newick string,depend on different species       """
        species_name = self.species
        fun_built_tree = getattr(SSTree, species_name)
        return fun_built_tree(genes_share_one_alignment)

    def write_batch_file(self, dot_input, dot_aln, hyphy_batch_file="", hyphy_result_file=""):
        """ write a .bf file for a alignment"""
        gene_id = os.path.basename(dot_input).split(os.path.extsep)[0]
        path_main = os.path.splitext(dot_input)[0]

        if "" == hyphy_batch_file:
            hyphy_batch_file = path_main + ".bf"

        if "" == hyphy_result_file:
            hyphy_result_file = path_main + ".result"

        if "" == self.batch_content:
            raise BFError

        # replace begins here
        batch_content, num_hits = re.subn(self.f_input, dot_input, self.batch_content)
        self._error_no_hit(num_hits)

        # partition is optional
        if (0, 0) == self.partition:
            if self.f_partition in batch_content:
                batch_content, num_hits = re.subn(self.f_partition, "", batch_content)
                self._error_no_hit(num_hits)

        else:
            batch_content, num_hits = re.subn(self.f_partition, "%d-%d" % self.partition, batch_content)
            self._error_no_hit(num_hits)

        batch_content, num_hits = re.subn(self.f_mdl, self.mdl_file, batch_content)
        self._error_no_hit(num_hits)

        # only support 1 matrix now :2014-5-26
        batch_content, num_hits = re.subn(self.f_matrix_name, self.matrix_name[0], batch_content)
        self._error_no_hit(num_hits)

        if self.use_given_tree:
            tree_newick_string = self.tree_definition_external
        else:
            genes_share_aln = pHdata.aln_reader(dot_aln)
            tree_newick_string = self.build_tree(genes_share_aln)

        batch_content, num_hits = re.subn(self.f_tree, tree_newick_string, batch_content)
        self._error_no_hit(num_hits)

        batch_content, num_hits = re.subn(self.f_output, hyphy_result_file, batch_content)
        self._error_no_hit(num_hits)

        self.check_whether_incomplete(batch_content)

        with open(name=hyphy_batch_file, mode="w") as bf_writer:
            bf_writer.write(batch_content)

    def check_whether_incomplete(self, batch_content):
        """ """
        flag_in_file = self.f_output in batch_content or \
                       self.f_tree in batch_content or \
                       self.f_mdl in batch_content or \
                       self.f_input in batch_content
        if flag_in_file:
            raise ParameterNotComplete

    @staticmethod
    def _error_no_hit(num_hit):
        if 0 == num_hit:
            raise ReplaceError


###########################################################
# : define exceptions for bf handler
class BFError(Exception):
    """    basic error type for bf_handler    """
    pass


class ParameterNotComplete(BFError):
    """    lack of some parameter. unable to make a functional batch file    """
    pass


class ReplaceError(BFError):
    """    error happens when replace text in template batch file    """
    pass


if __name__ == "__main__":
    pass