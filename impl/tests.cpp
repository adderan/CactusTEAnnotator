#include <math.h>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <iostream>
#include <climits>

#include "PairwiseDistances.h"
#include "sonLib.h"

using namespace std;

/*

void kmerDistanceTest() {
    int kmerLength = 3;
	int numSeqs = 3;
    char *a = (char*)malloc(1000*sizeof(char));
    char *b = (char*)malloc(1000*sizeof(char));
	char *c = (char*)malloc(1000*sizeof(char));

    strcpy(a, "AGAGCATGTGACTATGTGTATGTCGATGTGATGCTGATGCTAGCTAGCTAGCTGATCGATGC\0");
    strcpy(b, "AGAGCATGCGGCATATCGATCGTAGCACTATGTGTATGTCGATGTGATGCTGATGCTAGCTAGCTAGCTGATCGATGC\0");

	//reverse complement of a
	strcpy(c, "GCATCGATCAGCTAGCTAGCTAGCATCAGCATCACATCGACATACACATAGTCACATGCTCT");

    char **seqs = (char**)malloc(numSeqs*sizeof(char*));
    seqs[0] = a;
    seqs[1] = b;
	seqs[2] = c;

    int numHashes = 20000;

    double **dist = getDistances(seqs, numSeqs, kmerLength, numHashes);

    double exactDist = exactJaccardDistance(a, b, kmerLength);

    cerr << "Length = " << strlen(a) << endl;
    cerr << "dist = " << dist[1][0] << endl;
	cerr << "dist a c = " << dist[2][0] << endl;
    cerr << "exact dist = " << exactDist << endl;

	//a and c are reverse complements
	assert(dist[2][0] == -1.0);
}
*/


int main(int argc, char **argv) {
    //kmerDistanceTest();
}
