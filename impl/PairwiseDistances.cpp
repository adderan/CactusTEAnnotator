#include <vector>
#include <string>
#include <stack>
#include <iostream>
#include <getopt.h>
#include <set>
#include <algorithm>

#include "MurmurHash3.h"
#include <boost/math/distributions/binomial.hpp>

#include "PairwiseDistances.h"

extern "C" {
#include "sonLib.h"
#include "bioioC.h"
#include "commonC.h"

}


using namespace std;
using namespace boost::math;

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

double jaccardPValue(int sharedKmers, int totalKmers, int a_length, int b_length, int kmerLength) {
    if (sharedKmers == 0) {
        return 1.0;
    }
    int alphabet_size = 4;
    double kmerSpace = pow(alphabet_size, kmerLength);

    double px = 1.0 / (1.0 + kmerSpace/a_length);
    double py = 1.0 / (1.0 + kmerSpace/b_length);

    double r = px * py / (px + py  - px * py);
    //cerr << "shared kmers = " << sharedKmers << endl;
    //cerr << "total kmers = " << totalKmers << endl;
    //cerr << "r = " << r << endl;

    return cdf(complement(binomial(totalKmers, r), sharedKmers - 1));
}

vector<tuple<int, int, double> > getDistancesExact(char **sequences, int numSeqs, int kmerLength) {
    toUpperCase(sequences, numSeqs);
    vector<tuple<int, int, double> > pValues;
    for (int i = 0; i < numSeqs; i++) {
        for (int j = 0; j < i; j++) {
            double pValue = exactJaccardDistance(sequences[i], sequences[j], kmerLength);
            pValues.push_back(make_tuple(i, j, pValue));
        }
    }
    return pValues;
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


double minhashJaccard(vector<uint32_t> &a, vector<uint32_t> &b, int length_a, int length_b, int kmerLength) {
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
    return jaccardPValue(matches, total, length_a, length_b, kmerLength);
}


set<set<long> > buildClusters(vector<tuple<int, int, double> > &pValues, int nSeqs, double confidenceLevel) {

    //sort sequence pairs by ascending p-value
    sort(pValues.begin(), pValues.end(), [](tuple<int, int, double> a, tuple<int, int, double> b) {
        return (get<2>(a) < get<2>(b));
    });

    stUnionFind *components = stUnionFind_construct();
    for (int i = 0; i < nSeqs; i++) {
        stUnionFind_add(components, (void*)i+1);
    }
    for (auto &t: pValues) {
        int i = get<0>(t);
        int j = get<1>(t);
        double pValue = get<2>(t);

        int64_t size_i = stUnionFind_getSize(components, (void*)i+1);
        int64_t size_j = stUnionFind_getSize(components, (void*)j+1);

        double adjustedThreshold = confidenceLevel/(size_i * size_j);
        if (pValue < adjustedThreshold) {
            stUnionFind_union(components, (void*) i+1, (void*) j+1);
        }

    }

    set<set<long> > partitioning;
    stUnionFindIt *it = stUnionFind_getIterator(components);
    stSet *component;
    while((component = stUnionFindIt_getNext(it)) != NULL) {
        void *node;
        stSetIterator *nodeIt = stSet_getIterator(component);
        set<long> cluster;
        while((node = stSet_getNext(nodeIt)) 
                != NULL) {
            cluster.insert((long)node - 1);
        }
        partitioning.insert(cluster);
        stSet_destructIterator(nodeIt);
    }
    stUnionFind_destructIterator(it);
    stUnionFind_destruct(components);
    return partitioning;
}

vector<tuple<int, int, double> > getDistances(char **seqs, int numSeqs, int kmerLength, int numHashes) {

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

    vector<tuple<int, int, double> > pValues;
    for (int i = 0; i < numSeqs; i++) {
        for (int j = 0; j < i; j++) {
            //cerr << "seq = " << string((char*)sequences->list[i]) << endl;
            double p = minhashJaccard(sketches[i], sketches[j], strlen(seqs[i]), strlen(seqs[j]), kmerLength);

            pValues.push_back(make_tuple(i, j, p));
        }
    }


    return pValues;
}

