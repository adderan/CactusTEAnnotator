#include <vector>
#include <string>
#include <stack>
#include <iostream>
#include <tuple>
#include <getopt.h>
#include <set>

#include "MurmurHash3.h"

extern "C" {
#include "sonLib.h"
#include "bioioC.h"
#include "commonC.h"
}

#include "PairwiseDistances.h"

using namespace std;

int main(int argc, char **argv) {
    char *sequencesFilename = NULL;
    int kmerLength = 10;
    int numHashes = 200;
    bool exact = false;
    double confidenceLevel = 0.05;

    /*
     * Parse the options.
     */
    int i;
    while (1) {
        static struct option long_options[] = { 
            { "sequences", required_argument, 0, 'a' }, 
            { "kmerLength", required_argument, 0, 'b' }, 
            { "numHashes", required_argument, 0, 'c'},
            { "confidenceLevel", required_argument, 0, 'd'},
            { "exact", no_argument, 0, 'e'},
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d:e", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                sequencesFilename = stString_copy(optarg);
                break;
            case 'b':
                i = sscanf(optarg, "%d", &kmerLength);
                assert(i == 1);
                break;
            case 'c':
                i = sscanf(optarg, "%d", &numHashes);
                assert(i == 1);
                break;
            case 'd':
                i = sscanf(optarg, "%lf", &confidenceLevel);
                assert(i == 1);
                break;
            case 'e':
                exact = true;
                break;
            default:
                return 1;
        }
    }

    struct List *sequences = constructEmptyList(0, NULL);
    struct List *seqNames = constructEmptyList(0, free);
    struct List *seqLengths = constructEmptyList(0, free);
    FILE *sequencesFile = fopen(sequencesFilename, "r");
    fastaRead(sequencesFile, sequences, seqLengths, seqNames);
    cerr << "Read " << seqLengths->length << " sequences" << endl;

    vector<tuple<int, int, double> > pValues;
    if (exact) {
        pValues = getDistancesExact((char**)sequences->list, sequences->length, kmerLength);
    }
    else {
	    pValues = getDistances((char**)sequences->list, sequences->length, kmerLength, numHashes);
    }


    set<set<long> > partitioning = buildClusters(pValues, sequences->length, confidenceLevel);
    for (auto &cluster: partitioning) {
        for (auto &seqNum: cluster) {
            printf("%s ", seqNames->list[seqNum]);
        }
        printf("\n");
    }
}
