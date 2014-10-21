#include "CUSTOM_MODEL_FILE";

MAXIMUM_ITERATIONS_PER_VARIABLE = 10000000000;

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

function BuildCodonFrequencies (obsF)
{
	PIStop = 1.0;
	result = {61,1};
	hshift = 0;

	for (h=0; h<64; h=h+1)
	{
		first = h$16;
		second = h%16$4;
		third = h%4;
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
			PIStop = PIStop-obsF[first][0]*obsF[second][1]*obsF[third][2];
			continue;
		}
		result[h-hshift][0]=obsF[first][0]*obsF[second][1]*obsF[third][2];
	}
	return result*(1.0/PIStop);
}

DataSet nucleotideSequences = ReadDataFile ("INPUT_FILE");
DataSetFilter filteredData = CreateFilter (nucleotideSequences, 3, "PARTITION_DEFINITION", "", "TAA,TAG,TGA");
HarvestFrequencies (observedFreq,filteredData,3,1,1);
vectorOfFrequencies = BuildCodonFrequencies (observedFreq);
Model MyModel = (alphaWPsiMatrixP, vectorOfFrequencies);

Tree myTree = TREE_NEWICK_STRING;


LikelihoodFunction  theLnLik = (filteredData, myTree);
Optimize (paramValues, theLnLik);


fprintf ("OUTPUT_RESULT_FILE", CLEAR_FILE);
fprintf ("OUTPUT_RESULT_FILE", "Equilibrium codon frequency\n");
for (i=0;i<61;i=i+1)
{
	fprintf ("OUTPUT_RESULT_FILE", vectorOfFrequencies[i][0], "\n");
}
fprintf ("OUTPUT_RESULT_FILE", "Q matrix\n");
for (i=0;i<61;i=i+1)
{
	for (j=0;j<61;j=j+1)
	{
		fprintf ("OUTPUT_RESULT_FILE", alphaWPsiMatrixP[i][j], "\t");
	}
	fprintf ("OUTPUT_RESULT_FILE", "\n");
}
fprintf ("OUTPUT_RESULT_FILE", "Parameters\n");
fprintf ("OUTPUT_RESULT_FILE", "aAC = ", aAC, "\n");
fprintf ("OUTPUT_RESULT_FILE", "aAG = ", aAG, "\n");
fprintf ("OUTPUT_RESULT_FILE", "aAT = ", aAT, "\n");
fprintf ("OUTPUT_RESULT_FILE", "aCG = ", aCG, "\n");
fprintf ("OUTPUT_RESULT_FILE", "aCT = ", aCT, "\n");
fprintf ("OUTPUT_RESULT_FILE", "aGT = ", aGT, "\n");
fprintf ("OUTPUT_RESULT_FILE", "w = ", w, "\n");
fprintf ("OUTPUT_RESULT_FILE", "psi = ", psi, "\n");
fprintf ("OUTPUT_RESULT_FILE", theLnLik, "\n");



fprintf ("OUTPUT_RESULT_FILE", "\n number of local parmeters ---->	", paramValues[1][1]);
fprintf ("OUTPUT_RESULT_FILE", "\n number of global parmeters ---->	", paramValues[1][2]);

lnL = paramValues[1][0];
np = paramValues[1][1];
AIC = -2*lnL+2*np;

fprintf ("OUTPUT_RESULT_FILE", "\n AIC value ---->	", AIC);


