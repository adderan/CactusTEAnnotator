#include <vector>
#include <string>
#include <stack>
#include <iostream>
#include <tuple>
#include <getopt.h>
#include "MurmurHash3.h"

extern "C" {
#include "sonLib.h"
#include "bioioC.h"
#include "commonC.h"
}

#include "Minhash.h"

using namespace std;

char **readSequences(char *sequencesFilename, char ***sequenceNames, int *nSeqs) {
    struct List *sequences = constructEmptyList(0, NULL);
    struct List *seqNames = constructEmptyList(0, free);
    struct List *seqLengths = constructEmptyList(0, free);
    FILE *sequencesFile = fopen(sequencesFilename, "r");
    fastaRead(sequencesFile, sequences, seqLengths, seqNames);
    fclose(sequencesFile);
    *nSeqs = seqLengths->length;
    *sequenceNames = (char**)seqNames->list;
    return (char**)sequences->list;
}

int main(int argc, char **argv) {
    char *sequencesFilename = NULL;
    int kmerLength = 10;
    //int numHashes = 200;
    bool exact = false;
	//bool distancesOnly = false;
    char *clusterAFilename = NULL;
    char *clusterBFilename = NULL;
	int sketchSize = 100;

    /*
     * Parse the options.
     */
    while (1) {
        static struct option long_options[] = { 
            { "sequences", required_argument, 0, 'a' }, 
            { "kmerLength", required_argument, 0, 'b' }, 
            { "exact", no_argument, 0, 'd'},
            { "clusterA", required_argument, 0, 'f'},
            { "clusterB", required_argument, 0, 'g'},
			{ "sketchSize", required_argument, 0, 'h'},
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:d:f:g:h:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                sequencesFilename = stString_copy(optarg);
                break;
            case 'b':
                sscanf(optarg, "%d", &kmerLength);
                break;
            case 'd':
                exact = true;
                break;
            case 'f':
                clusterAFilename = stString_copy(optarg);
                break;
            case 'g':
                clusterBFilename = stString_copy(optarg);
                break;
			case 'h':
				sscanf(optarg, "%d", &sketchSize);
				break;
            default:
                return 1;
        }
    }

    if (sequencesFilename) {
        int numSequences;
        char **seqNames;
        char **sequences = readSequences(sequencesFilename, &seqNames, &numSequences);
        cerr << "Read " << numSequences << " sequences" << endl;

        vector<tuple<int, int, double> > distances;
        if (exact) {
            distances = getDistancesExact(sequences, numSequences, kmerLength);
        }
        else {
            distances = getDistances(sequences, numSequences, kmerLength, sketchSize);
        }

		for (auto &t: distances) {
			printf("%s %s %f\n", seqNames[get<0>(t)], seqNames[get<1>(t)], get<2>(t));
		}
    }

    else if (clusterAFilename && clusterBFilename) {
        char **clusterASeqNames;
        int clusterANumSeqs;
        char **clusterASequences = readSequences(clusterAFilename, &clusterASeqNames, &clusterANumSeqs);

        char **clusterBSeqNames;
        int clusterBNumSeqs;
        char **clusterBSequences = readSequences(clusterBFilename, &clusterBSeqNames, &clusterBNumSeqs);

        double distance = getDistanceBetweenFamilies(clusterASequences, clusterBSequences, clusterANumSeqs, clusterBNumSeqs, kmerLength, sketchSize);

        cout << distance << endl;

    }

}
