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

double exactJaccardDistance(char *a, char*b, int kmerLength) {
    stSet *sharedKmers = stSet_construct();
    stSet *totalKmers = stSet_construct();

    for (int i = 0; i < strlen(a) - kmerLength; i++) {
        uint32_t hash = hashKmer(a + i, kmerLength);
        stSet_insert(totalKmers, (void*)hash);
    }
    for (int i = 0; i < strlen(b) - kmerLength; i++) {
        uint32_t hash = hashKmer(b + i, kmerLength);
        stSet_insert(totalKmers, (void*)hash);
    }

    for (int i = 0; i < strlen(a) - kmerLength; i++) {
        for (int j = 0; j < strlen(b) - kmerLength; j++) {
            if (strncmp(a + i, b + j, kmerLength) == 0) {
                uint32_t hash = hashKmer(a + i, kmerLength);
                stSet_insert(sharedKmers, (void*)hash);
                break;
            }
        }
    }
    cerr << "Matching kmers = " << stSet_size(sharedKmers) << endl;
    cerr << "Total kmers = " << stSet_size(totalKmers) << endl;
    return stSet_size(sharedKmers)/(double)stSet_size(totalKmers);
}

void kmerDistanceTest() {
    int kmerLength = 3;
    char *a = (char*)malloc(1000*sizeof(char));
    char *b = (char*)malloc(1000*sizeof(char));

    strcpy(a, "AGAGCATGTGACTATGTGTATGTCGATGTGATGCTGATGCTAGCTAGCTAGCTGATCGATGC\0");
    strcpy(b, "AGAGCATGCGGCATATCGATCGTAGCACTATGTGTATGTCGATGTGATGCTGATGCTAGCTAGCTAGCTGATCGATGC\0");

    char **seqs = (char**)malloc(2*sizeof(char*));
    seqs[0] = a;
    seqs[1] = b;

    int numHashes = 20000;

    uint32_t **minhashValues = precompute_minhash(seqs, 2, kmerLength, numHashes);

    double dist = minhash_jaccard(0, 1, minhashValues, numHashes);

    double exactDist = exactJaccardDistance(a, b, kmerLength);
    cerr << "Length = " << strlen(a) << endl;
    cerr << "dist = " << dist << endl;
    cerr << "exact dist = " << exactDist << endl;

}


int main(int argc, char **argv) {
    kmerDistanceTest();
}
