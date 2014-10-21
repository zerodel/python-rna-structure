
_Genetic_Code = {

		{14,/*AAA*/ 13,/*AAC*/ 14,/*AAG*/  13,/*AAT*/
		  7, /*ACA*/ 7, /*ACC*/ 7, /*ACG*/  7, /*ACT*/
		 19, /*AGA*/ 5, /*AGC*/19, /*AGG*/  5, /*AGT*/
		  2, /*ATA*/ 2, /*ATC*/	3, /*ATG*/  2, /*ATT*/
		 12,/*CAA*/ 11,/*CAC*/ 12,/*CAG*/  11,/*CAT*/
		  6, /*CCA*/ 6, /*CCC*/ 6, /*CCG*/  6, /*CCT*/
		 19,/*CGA*/ 19,/*CGC*/ 19,/*CGG*/  19,/*CGT*/
		  1, /*CTA*/ 1, /*CTG*/ 1, /*CTC*/  1, /*CTT*/
		 16,/*GAA*/ 15,/*GAC*/ 16,/*GAG*/  15,/*GAT*/
		  8, /*GCA*/ 8, /*GCC*/ 8, /*GCG*/  8, /*GCT*/
		 20,/*GGA*/ 20,/*GGC*/ 20,/*GGG*/  20,/*GGT*/
		  4, /*GTA*/ 4, /*GTC*/ 4, /*GTG*/  4, /*GTT*/
		 10,/*TAA*/  9, /*TAC*/10,/*TAG*/   9, /*TAT*/
		  5, /*TCA*/ 5, /*TCC*/ 5, /*TCG*/  5, /*TCT*/
		 10,/*TGA*/ 17,/*TGC*/ 18,/*TGG*/  17,/*TGT*/
		  1, /*TTA*/ 0, /*TTC*/ 1, /*TTG*/  0  /*TTT*/ }
	};

GeneticCodeExclusions = "TAA,TAG,TGA";

DataSet nucleotideSequences = ReadDataFile ("INPUT_FILE");

fprintf (stdout,"\n ======data set is \n", ds);
fprintf ("OUTPUT_RESULT_FILE", CLEAR_FILE);

DataSetFilter filteredData = CreateFilter (nucleotideSequences,3,"PARTITION_DEFINITION","", GeneticCodeExclusions);

fprintf (stdout, "\n now filter is ok");
fprintf ("OUTPUT_RESULT_FILE",filteredData);

HarvestFrequencies (FEV, filteredData, 1, 1, 0);
fprintf (stdout, "\n frequency ready \n");

ACCEPT_ROOTED_TREES = 1;
ModelGeneticCode = _Genetic_Code;

modelType := 0;
SKIP_MODEL_PARAMETER_LIST := 1;

#include "CUSTOM_MODEL_FILE";

Tree givenTree = TREE_NEWICK_STRING;

fprintf(stdout , "\n model included\n define model now ......\n");
fprintf ("OUTPUT_RESULT_FILE", "\n =====Newick String===== \n",givenTreen);
fprintf (stdout, "\n tree defined \n");

LikelihoodFunction lf = (filteredData,givenTree);

fprintf (stdout, "\n likelihood function is ok \n");

Optimize	(res,lf);

LIKELIHOOD_FUNCTION_OUTPUT = 1;
fprintf (stdout, "\n______________RESULTS______________\n",lf);
fprintf ("OUTPUT_RESULT_FILE","\n______________RESULTS______________\n",lf);
lnL = res[1][0];
np = res[1][1];

fprintf (stdout, "\n likelihood is ----> ", lnL);
fprintf (stdout, "\n number of local parmeters ---->	", res[1][1]);
fprintf (stdout, "\n number of global parmeters ---->	", res[1][2]);
fprintf ("OUTPUT_RESULT_FILE", "\n number of local parmeters ---->	", res[1][1]);
fprintf ("OUTPUT_RESULT_FILE", "\n number of global parmeters ---->	", res[1][2]);

AIC = -2*lnL+2*np;

fprintf (stdout, "\n AIC value ---->	",AIC);
fprintf ("OUTPUT_RESULT_FILE", "\n AIC value ---->	", AIC);

