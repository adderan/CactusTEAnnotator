#include <vector>
#include <string>
#include <stack>
#include <iostream>
#include <getopt.h>

#include "sonLib.h"
#include "bioioC.h"
#include "commonC.h"
#include "MurmurHash3.h"

#include "PairwiseDistances.h"


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

char *reverseComplement(char *seq) {
	char *reverse = (char*) malloc(sizeof(char)*strlen(seq));
	for (int i = 0; i < strlen(seq); i++) {
		char base = seq[strlen(seq) - i - 1];
		if (base == 'A') {
			reverse[i] = 'T';
		}
		else if (base == 'T') {
			reverse[i] = 'A';
		}
		else if (base == 'G') {
			reverse[i] = 'C';
		}
		else if (base == 'C') {
			reverse[i] = 'G';
		}

	}
	return reverse;
}

uint32_t **precompute_minhash(char **seqs, int numSeqs, int kmerLength, int numHashes, uint32_t *a, uint32_t *b, uint32_t p) {
    uint32_t **minhashValues = (uint32_t**)malloc(sizeof(uint32_t*)*numSeqs);


    for (int j = 0; j < numSeqs; j++) {

        minhashValues[j] = (uint32_t*) malloc(sizeof(uint32_t)*numHashes);

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
            minhashValues[j][i] = min;
        }
        free(hashStart);
    }

    return minhashValues;
}

double minhashJaccard(uint32_t *values_a, uint32_t *values_b, int numHashes) {
    int matches = 0;
    for (int i = 0; i < numHashes; i++) {
        if (values_a[i] == values_b[i]) {
            matches++;
        }
    }
    return (double) matches/(double) numHashes;
}

double **getDistances(char **seqs, int numSeqs, int kmerLength, int numHashes) {

    //convert to uppercase
    for (int i = 0; i < numSeqs; i++) {
        char *seq = seqs[i];
        for (int k = 0; k < strlen(seq); k++) {
            seq[k] = (char)toupper(seq[k]);
        }
    }

	char **revSeqs = (char**)malloc(sizeof(char*)*numSeqs);
	for (int i = 0; i < numSeqs; i++) {
		revSeqs[i] = reverseComplement(seqs[i]);
	}

    uint32_t p = (1 << 16) - 1;
    uint32_t *a = (uint32_t*) malloc(sizeof(uint32_t)*numHashes);
    uint32_t *b = (uint32_t*) malloc(sizeof(uint32_t)*numHashes);
    for (int i = 0; i < numHashes; i++) {
        a[i] = rand() * (float) p;
        b[i] = rand() * (float) p;
    }

    uint32_t **minhashValues = 
		precompute_minhash(seqs, numSeqs, kmerLength, numHashes, a, b, p);

	uint32_t **revMinhashValues = 
		precompute_minhash(revSeqs, numSeqs, kmerLength, numHashes, a, b, p);

	double **distances = (double**) malloc(sizeof(double*) * numSeqs);
    for (int i = 0; i < numSeqs; i++) {
		distances[i] = (double*) malloc(sizeof(double) * numSeqs);
        for (int j = 0; j < i; j++) {
            //cerr << "seq = " << string((char*)sequences->list[i]) << endl;
            double dist = minhashJaccard(minhashValues[i], minhashValues[j], numHashes);
			double revDist = minhashJaccard(minhashValues[i], revMinhashValues[j], numHashes);

			if (dist > revDist) {
				distances[i][j] = dist;
			}
			else {
				distances[i][j] = -1*revDist;
			}
        }
    }
	return distances;
}

