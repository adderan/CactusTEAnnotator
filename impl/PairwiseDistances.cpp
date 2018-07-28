#include <vector>
#include <string>
#include <stack>
#include <iostream>
#include <getopt.h>

#include "sonLib.h"
#include "bioioC.h"
#include "commonC.h"
#include "MurmurHash3.h"


using namespace std;
//using namespace boost::numeric::ublas;

uint32_t hashKmer(char *seq, int length, int seed) {
    for (int i = 0; i < length; i++) {
        if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'T' && seq[i] != 'G') return -1;
    }
    char data[8];
    MurmurHash3_x86_32(seq, length, seed, data);
    return *(uint32_t*)data;
}

uint32_t minhash(char *seq, int kmerLength, int seed) {
    uint32_t min = INT_MAX;
    //cerr << "seq = " << string(seq) << endl;
    for (int i = 0; i < strlen(seq) - kmerLength; i++) {
        uint32_t hash = hashKmer(seq + i, kmerLength, seed);
        //cerr << "hash = " << hash << endl;
        if (hash < min) {
            min = hash;
        }
    }
    return min;
}

uint32_t **precompute_minhash(char **seqs, int numSeqs, int kmerLength, uint32_t *seeds, int numSeeds) {
    uint32_t **minhashValues = (uint32_t**)malloc(sizeof(uint32_t*)*numSeeds);
    for (int i = 0; i < numSeeds; i++) {
        minhashValues[i] = (uint32_t*) malloc(sizeof(uint32_t)*numSeqs);
        for (int j = 0; j < numSeqs; j++) {
            minhashValues[i][j] = minhash((char*)seqs[j], kmerLength, seeds[i]);
        }
    }
    return minhashValues;
}


double minhash_jaccard(int a, int b, uint32_t **minhashValues, int numSeeds) {
    int matches = 0;
    for (int i = 0; i < numSeeds; i++) {
        if (minhashValues[i][a] == minhashValues[i][b]) {
            matches++;
        }
    }
    return (double) matches/(double) numSeeds;
}
