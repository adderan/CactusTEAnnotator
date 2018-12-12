#include <math.h>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <iostream>
#include <climits>

#include "Minhash.h"
#include "sonLib.h"

using namespace std;

void kmerDistanceTest() {
    int kmerLength = 5;
	int numSeqs = 3;
    char *a = (char*)malloc(1000*sizeof(char));
    char *b = (char*)malloc(1000*sizeof(char));
	char *c = (char*)malloc(1000*sizeof(char));

    strcpy(a, "AGAGCATGTGACTATGTGTATGTCGATGTGATGCTGATGCTAGCTAGCTAGCTGATCGATGC\0");
    strcpy(b, "AGAGCATGCGGCATATCGATCGTAGCACTATGTGTATGTCGATGTGATGCTGATGCTAGCTAGCTAGCTGATCGATGC\0");

	//reverse complement of a
	strcpy(c, "GCATCGATCAGCTAGCTAGCTAGCATCAGCATCACATCGACATACACATAGTCACATGCTCT");


    vector<uint32_t> a_sketch = buildSketch(a, kmerLength);
    vector<uint32_t> b_sketch = buildSketch(b, kmerLength);

    double minhashDist = minhashJaccard(a_sketch, b_sketch, strlen(a), strlen(b), kmerLength);

    double exactDist = exactJaccardDistance(a, b, kmerLength);

    cerr << "Length = " << strlen(a) << endl;
    cerr << "dist a b = " << minhashDist << endl;
    cerr << "exact dist = " << exactDist << endl;

}


int main(int argc, char **argv) {
    kmerDistanceTest();
}
