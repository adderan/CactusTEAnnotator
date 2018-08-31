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

int main(int argc, char **argv) {
    char *sequencesFilename = NULL;
    int kmerLength = 10;
    int numHashes = 200;

    /*
     * Parse the options.
     */
    int i;
    while (1) {
        static struct option long_options[] = { 
            { "sequences", required_argument, 0, 'a' }, 
            { "kmerLength", required_argument, 0, 'b' }, 
            { "numHashes", required_argument, 0, 'c'},
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:", long_options, &option_index);

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

	double **distances = getDistances((char**)sequences->list, sequences->length, kmerLength, numHashes);

	for (int i = 0; i < sequences->length; i++) {
		for (int j = 0; j < i; j++) {

			if (distances[i][j] > 0.0) {
				printf("%s %s %f %c\n", seqNames->list[i], seqNames->list[j], distances[i][j], '+');
			}

			else {
				printf("%s %s %f %c\n", seqNames->list[i], seqNames->list[j], -1.0 * distances[i][j], '-');
			}

		}
	}
}
