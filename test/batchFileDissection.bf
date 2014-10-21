global aAC = 1;
global aAG = 1;
global aAT = 1;
global aCG = 1;
global aCT = 1;
global aGT = 1;
global w =1;
global psi =1;

ACCEPT_ROOTED_TREES = 1;
LIKELIHOOD_FUNCTION_OUTPUT = 1;

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
codon_list_hypothesis =  {
		      {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA",
                         "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC",
                         "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG",
                         "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC",
                         "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"}
			 /*{AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT, AGA, AGC, AGG, AGT, ATA, ATC, ATG, ATT, CAA, CAC, CAG, CAT, CCA, CCC, CCG, CCT, CGA, CGC, CGG, CGT, CTA, CTC, CTG, CTT, GAA, GAC, GAG, GAT, GCA, GCC, GCG, GCT, GGA, GGC, GGG, GGT, GTA, GTC, GTG, GTT, TAC, TAT, TCA, TCC, TCG, TCT, TGC, TGG, TGT, TTA, TTC, TTG, TTT }*/

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

DataSet nucleotideSequences = ReadDataFile ("tmp.input");
DataSetFilter filteredData = CreateFilter (nucleotideSequences, 3, "", "", "TAA,TAG,TGA");
HarvestFrequencies (observedFreq,filteredData,3,1,1);
vectorOfFrequencies = BuildCodonFrequencies (observedFreq);

fprintf (stdout,"\n===============================>raw \n", nucleotideSequences,"\n");
fprintf (stdout, "\n==============================>filtered raw ,\n", filteredData,"\n");
fprintf (stdout, "\n==============================>obf \n", observedFreq, "\n");

/*
for (i=0;i<61;i=i+1)
{
	fprintf (stdout, codon_list_hypothesis[0][i], "\t", vectorOfFrequencies[i][0], "\n");
}
*/


/* define GTR model  */
#include "rebuild_model.mdl";
Model myModel = (alphaWPsiMatrixP, vectorOfFrequencies);

Tree myTree = (g1,g2);


/*
Tree myTree = (((macaque:0.02602176453689536428,(chimp:0.00432283016370958121,human:0.00384313582540942496):0.00915641286880421179):0.03904538064229965549,(dog:0.06293344444628010126,(cow:0.07356214799052007702,horse:0.04838456774727731280):0.00508183725708111558):0.01752218773337325258):0.09364969675060515197,rat:0.03601197114550155204,mouse:0.02852390293940118213):0.0;
*/

LikelihoodFunction  theLnLik = (filteredData, myTree);


/*wait here, */

/*
Optimize (res, theLnLik);


fprintf ("ENSMUSP00000029836gtr.result", CLEAR_FILE);
fprintf ("ENSMUSP00000029836gtr.result", "Equilibrium codon frequency\n");
for (i=0;i<61;i=i+1)
{
	fprintf ("ENSMUSP00000029836gtr.result", vectorOfFrequencies[i][0], "\n");
}
fprintf ("ENSMUSP00000029836gtr.result", "Q matrix\n");
for (i=0;i<61;i=i+1)
{
	for (j=0;j<61;j=j+1)
	{
		fprintf ("ENSMUSP00000029836gtr.result", alphaWPsiMatrixP[i][j], "\t");
	}
	fprintf ("ENSMUSP00000029836gtr.result", "\n");
}

fprintf ("ENSMUSP00000029836gtr.result", theLnLik, "\n");

lnL = res[1][0];
np = res[1][1];

fprintf (stdout, "\n likelihood is ----> ", lnL);
fprintf (stdout, "\n number of local parmeters ---->	", res[1][1]);
fprintf (stdout, "\n number of global parmeters ---->	", res[1][2]);
fprintf ("ENSMUSP00000029836gtr.result", "\n number of local parmeters ---->	", res[1][1]);
fprintf ("ENSMUSP00000029836gtr.result", "\n number of global parmeters ---->	", res[1][2]);

AIC = -2*lnL+2*np;

fprintf (stdout, "\n AIC value ---->	",AIC);
fprintf ("ENSMUSP00000029836gtr.result", "\n AIC value ---->	", AIC);
*/