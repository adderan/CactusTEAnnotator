#include <vector>
#include <string>
#include <stack>
#include <iostream>
#include <getopt.h>
#include <set>
#include <algorithm>

#include "sonLib.h"
#include "bioioC.h"
#include "commonC.h"
#include "MurmurHash3.h"

#include "PairwiseDistances.h"


using namespace std;

void toUpperCase(char **seqs, int numSeqs) {
    for (int i = 0; i < numSeqs; i++) {
        char *seq = seqs[i];
        for (int k = 0; k < strlen(seq); k++) {
            seq[k] = (char)toupper(seq[k]);
        }
    }
}

uint32_t hashKmer(char *seq, int length) {
    for (int i = 0; i < length; i++) {
        if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'T' && seq[i] != 'G') return -1;
    }
    char data[8];
    MurmurHash3_x86_32(seq, length, 0, data);
    return *(uint32_t*)data;
}


double **getDistancesExact(char **sequences, int numSeqs, int kmerLength) {
    toUpperCase(sequences, numSeqs);
    double **distances = (double**) malloc(sizeof(double*) * numSeqs);
    for (int i = 0; i < numSeqs; i++) {
        distances[i] = (double*) malloc(sizeof(double) * numSeqs);
        for (int j = 0; j < i; j++) {
            distances[i][j] = exactJaccardDistance(sequences[i], sequences[j], kmerLength);
        }
    }
    return distances;
}

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

    int nSharedKmers = 0;
    for (int i = 0; i < strlen(a) - kmerLength; i++) {
        for (int j = 0; j < strlen(b) - kmerLength; j++) {
            if (strncmp(a + i, b + j, kmerLength) == 0) {
                /*
                fprintf(stderr, "A = ");
                for (int k = 0; k < kmerLength; k++) {
                    fprintf(stderr, "%c", a[i+k]);
                }
                fprintf(stderr, "\n");
                fprintf(stderr, "B = ");
                for (int k = 0; k < kmerLength; k++) {
                    fprintf(stderr, "%c", b[j+k]);
                }
                fprintf(stderr, "\n\n");
                */
                nSharedKmers += 2;
                uint32_t hash = hashKmer(a + i, kmerLength);
                stSet_insert(sharedKmers, (void*)hash);
                break;
            }
        }
    }
    fprintf(stderr, "distance = %f\n", nSharedKmers/((double)(strlen(a) - kmerLength) + (double)(strlen(b) - kmerLength)));
    fprintf(stderr, "Shared kmers %d\n", stSet_size(sharedKmers));
    fprintf(stderr, "Total kmers %d\n", stSet_size(totalKmers));
    return stSet_size(sharedKmers)/(double)stSet_size(totalKmers);
}


double minhashJaccard(vector<uint32_t> &a, vector<uint32_t> &b) {
    int matches = 0;
    int total = 0;
    int i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a.at(i) == b.at(j)) {
            i++;
            j++;
            matches++;
        }
        else if (a.at(i) < b.at(j)) {
            i++;
        }
        else if (a.at(i) > b.at(j)) {
            j++;
        }
        total++;
    }

    return (double) matches/(double) total;
}

double **getDistances(char **seqs, int numSeqs, int kmerLength, int numHashes) {

    toUpperCase(seqs, numSeqs);

    //Build minhash sketches
    vector<vector<uint32_t> > sketches(numSeqs);
    for (int i = 0; i < numSeqs; i++) {
        for (int j = 0; j < strlen(seqs[i]) - kmerLength; j++) {
            uint32_t hash = hashKmer(seqs[i] + j, kmerLength);
            sketches[i].push_back(hash);
        }
        sort(sketches[i].begin(), sketches[i].end());

    }

	double **distances = (double**) malloc(sizeof(double*) * numSeqs);
    for (int i = 0; i < numSeqs; i++) {
		distances[i] = (double*) malloc(sizeof(double) * numSeqs);
        for (int j = 0; j < i; j++) {
            //cerr << "seq = " << string((char*)sequences->list[i]) << endl;
            double dist = minhashJaccard(sketches[i], sketches[j]);

			distances[i][j] = dist;
        }
    }
	return distances;
}

