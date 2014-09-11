__author__ = 'root'

import unittest

import pyRNAsnp.pyRNAsnp as RNAsnp




class MyTestCase(unittest.TestCase):
    def test_something(self):
        print "hello world"
        self.assertEqual(True, True)


    def test_get_all_from_one_file(self):
        filename_rnasnp = "rnasnp4test.rnasnp"
        rnasnp_content,rnasnp_head = RNAsnp.rnasnp_extractor(filename_rnasnp)
        self.assertEqual(["A1G"], RNAsnp.get_snp_list(rnasnp_content))

    def test_mutation_site(self):
        mutation1 = "A"
        pos = 2
        snps_in_mind = ["A2G", "A2C", "A2T"]
        snp_output = RNAsnp.mutation_one_site(mutation1, pos)
        print str(snp_output)
        for snp in snps_in_mind:
            self.assertIn(snp, snp_output)

    def test_rebuild_sequence(self):
        pass
        #give a sequence
        seq_given = "ATGTGTGCCCTTGAAACCCC"
        # here we build fake snp file
        tmp_file_name = "./fake_snp.rnasnp"
        tmp_file_with_line = "./snp_line.rnasnp"
        table_head = "SNP	w	Slen	GC	interval	d	pvalue1	ewin	interval	d_max	pvalue2"
        describe_line_one = "200	2463	0.5373	1-26	0.0361	0.0943	200	1-53	0.0072	0.4643"
        snp_line_picked_out = "200	2463	0.5373	1-26	0.0361	0.0943	200	1-53	0.0072	0.0643"
        p2_less = 0.0643
        line_num_pick = 9
        with open(tmp_file_name, "w") as snp_faker:
            snp_faker.write(table_head)

            snps = []
            for index_nt, nt in enumerate(seq_given):
                snps.extend(RNAsnp.mutation_one_site(nt, index_nt+1))

            snp_faker.write("\n".join(snps))

        with open(tmp_file_with_line, "w") as snps_add_line:
            snps_add_line.write(table_head)

            snp_lines = []
            for snp in snps:
                if snps.index(snp) == line_num_pick:
                    snp_lines.append("\t".join([snp, snp_line_picked_out]))
                else:
                    snp_lines.append("\t".join([snp, describe_line_one]))

            snps_add_line.writelines("\n".join(snp_lines))

        # use this file as input for pyRNAsnp.get_snp_list.
        rebuilded_seq = RNAsnp.rebuild_seq(tmp_file_name)
        self.assertEqual(seq_given, rebuilded_seq)

        # test how to pick out a line with less p-value
        picked_collected = []
        with open(tmp_file_with_line, "r") as reader:
            for line in reader.readlines():
                parts = line.split()
                p2 = float(parts[-1])
                if p2 < 0.1:
                    print line
                    picked_collected.append(p2)

        self.assertIn(p2_less, picked_collected)



    def test_rnasnp_built_cpd(self):
        """
        check the performance after a filter added
        :return:
        """
        snp_file_name = "./snp_line.rnasnp"
        rnasnp_list, rnasnp_list_head = RNAsnp.rnasnp_extractor(snp_file_name)
        print len(rnasnp_list)
        seq_rebuilt = RNAsnp.rebuild_seq(snp_file_name)
        print seq_rebuilt
        # rnasnp_built_cpd(list_rnasnp, head_list_rnasnp, sequence_string, list_codon=gc.codon_list):
        codon_pair = RNAsnp.rnasnp_built_cpd(rnasnp_list, rnasnp_list_head, seq_rebuilt)

        for key1 in codon_pair.keys():
            if codon_pair[key1]:
                print key1, str(codon_pair[key1])



if __name__ == '__main__':
    unittest.main()
