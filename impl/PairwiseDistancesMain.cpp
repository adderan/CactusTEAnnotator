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
    double distanceThreshold = 0.1;
	bool distancesOnly = false;

    /*
     * Parse the options.
     */
    int i;
    while (1) {
        static struct option long_options[] = { 
            { "sequences", required_argument, 0, 'a' }, 
            { "kmerLength", required_argument, 0, 'b' }, 
            { "distanceThreshold", required_argument, 0, 'c'},
            { "exact", no_argument, 0, 'd'},
			{ "distancesOnly", no_argument, 0, 'e'},
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d", long_options, &option_index);

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
                i = sscanf(optarg, "%lf", &distanceThreshold);
                assert(i == 1);
                break;
            case 'd':
                exact = true;
                break;
			case 'e':
				distancesOnly = true;
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

    vector<tuple<int, int, double> > distances;
    if (exact) {
        distances = getDistancesExact((char**)sequences->list, sequences->length, kmerLength);
    }
    else {
	    distances = getDistances((char**)sequences->list, sequences->length, kmerLength);
    }

	if (distancesOnly) {
		for (auto &t: distances) {
			printf("%s %s %f\n", seqNames->list[get<0>(t)], seqNames->list[get<1>(t)], get<2>(t));
		}
		exit(0);
	}
    set<set<long> > partitioning = buildClusters(distances, sequences->length, distanceThreshold);
    for (auto &cluster: partitioning) {
        for (auto &seqNum: cluster) {
            printf("%s ", seqNames->list[seqNum]);
        }
        printf("\n\n");
    }
    cerr << "Found " << partitioning.size() << " clusters" << endl;
}
