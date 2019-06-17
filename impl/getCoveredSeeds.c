#include "sonLib.h"
#include <getopt.h>
#include <stdlib.h>

#include "sonLib.h"
#include "bioioC.h"
#include "pairwiseAlignment.h"

int main(int argc, char **argv) {
	char *sequencesFilename = NULL;
	char *alignmentsFilename = NULL;
	char *seedsFilename = NULL;
	while (1) {
        static struct option long_options[] = {
            { "sequences", required_argument, 0, 'a' }, 
			{ "alignments", required_argument, 0, 'b'},
			{ "seeds", required_argument, 0, 'c'},
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
				sequencesFilename = strdup(optarg);
                break;
            case 'b':
				alignmentsFilename = strdup(optarg);
                break;
			case 'c':
				seedsFilename = strdup(optarg);
				break;
            default:
                return 1;
        }
    }

	FILE *alignmentsFile = fopen(alignmentsFilename, "r");
	FILE *sequencesFile = fopen(sequencesFilename, "r");
	FILE *seedsFile = fopen(seedsFilename, "r");


	//Just keep track of the hash lastz assigns
	//to each seed instead of computing it from scratch
	stHash *seedToLastzHash = stHash_construct();
	stHash *seedCounts = stHash_construct();
	uint64_t lastzHash;
	char seed[50];
	uint64_t count;
	char *line = NULL;
	size_t nbytes;
	while(getline(&line, &nbytes, seedsFile) != EOF) {
		int ret = sscanf(line, "%lx/%[^:]: %lu", &lastzHash, seed, &count);
		if (ret != 3) continue;

		uint64_t *lastzHashPtr = malloc(sizeof(uint64_t));
		memcpy(lastzHashPtr, &lastzHash, sizeof(uint64_t));
		stHash_insert(seedToLastzHash, (void*) stHash_stringKey(seed), (void*) lastzHashPtr);

		uint64_t *countPtr = malloc(sizeof(uint64_t));
		memcpy(countPtr, &count, sizeof(uint64_t));
		stHash_insert(seedCounts, (void*) stHash_stringKey(seed), (void*) countPtr);
	}

	stHash *seqNameToSeq = stHash_construct();
	struct List *seqs = constructEmptyList(0, NULL);
    struct List *seqLengths = constructEmptyList(0, free);
    struct List *headers = constructEmptyList(0, free);
	fastaRead(sequencesFile, seqs, seqLengths, headers);
	for (int i = 0; i < seqs->length; i++) {
		stHash_insert(seqNameToSeq, (void*) stHash_stringKey(headers->list[i]), (void*) seqs->list[i]);
	}

	//Detect the seed length and ignored positions
	int seedLength = strlen(seed);
	char *seedPattern = seed;

	char *seedWindow = calloc(sizeof(char), seedLength + 1);
	struct PairwiseAlignment *alignment = NULL;
	while ((alignment = cigarRead(alignmentsFile)) != NULL) {
		char *seq = stHash_search(seqNameToSeq, (void*) stHash_stringKey(alignment->contig1));
		for (int i = 0; i < strlen(seq) - seedLength; i++) {
			memcpy(seedWindow, seq + i, seedLength);
			for (int pos = 0; pos < seedLength; pos++) {
				if (seedPattern[pos] == 'x') seedWindow[pos] = 'x';
			}
			uint64_t *seedWindowHash = (uint64_t*) stHash_search(seedToLastzHash, (void*) stHash_stringKey(seedWindow));
			if (!seedWindowHash) continue;
			uint64_t *seedMultiplicity = (uint64_t*) stHash_search(seedCounts, (void*) stHash_stringKey(seedWindow));

			printf("%.6lx\t%lu\t%lu\n", *seedWindowHash, *seedMultiplicity, strlen(seq) - seedLength);
		}
	}

	fclose(alignmentsFile);
	fclose(seedsFile);
	fclose(sequencesFile);
}
