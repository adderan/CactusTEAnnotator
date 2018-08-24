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

uint32_t hashKmer(char *seq, int length) {
    for (int i = 0; i < length; i++) {
        if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'T' && seq[i] != 'G') return -1;
    }
    char data[8];
    MurmurHash3_x86_32(seq, length, 0, data);
    return *(uint32_t*)data;
}

uint32_t **precompute_minhash(char **seqs, int numSeqs, int kmerLength, int numHashes) {
    uint32_t **minhashValues = (uint32_t**)malloc(sizeof(uint32_t*)*numHashes);

    uint32_t p = (1 << 16) - 1;
    uint32_t *a = (uint32_t*) malloc(sizeof(uint32_t)*numHashes);
    uint32_t *b = (uint32_t*) malloc(sizeof(uint32_t)*numHashes);
    for (int i = 0; i < numHashes; i++) {
        minhashValues[i] = (uint32_t*) malloc(sizeof(uint32_t)*numSeqs);
        a[i] = rand() * (float) p;
        b[i] = rand() * (float) p;
    }

    for (int j = 0; j < numSeqs; j++) {

        uint32_t *hashStart = (uint32_t*) malloc(sizeof(uint32_t)*(strlen(seqs[j]) - kmerLength));
        for (int k = 0; k < strlen(seqs[j]) - kmerLength; k++) {
            hashStart[k] = hashKmer(seqs[j] + k, kmerLength);
        }
        for (int i = 0; i < numHashes; i++) {

            uint32_t min = INT_MAX;
            int min_k = 0;
            for (int k = 0; k < strlen(seqs[j]) - kmerLength; k++) {

                uint32_t hash = (a[i]*hashStart[k] + b[i]) % p;
                
                if (hash < min) {
                    min = hash;
                    min_k = k;
                }
            }
            minhashValues[i][j] = min;
        }
        free(hashStart);
    }

    free(a);
    free(b);
    return minhashValues;
}


double minhash_jaccard(int a, int b, uint32_t **minhashValues, int numHashes) {
    int matches = 0;
    for (int i = 0; i < numHashes; i++) {
        if (minhashValues[i][a] == minhashValues[i][b]) {
            matches++;
        }
    }
    return (double) matches/(double) numHashes;
}
