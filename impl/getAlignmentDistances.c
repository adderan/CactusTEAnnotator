#include <stdio.h>
#include <getopt.h>

#include "bioioC.h"
#include "pairwiseAlignment.h"

int main(int argc, char **argv) {
	char *sequencesFilename = NULL;
	char *alignmentsFilename = NULL;
	while (1) {
		static struct option long_options[] = {
			{ "alignments", required_argument, 0, 'a' }, 
            { "sequences", required_argument, 0, 'b'},
			{ 0, 0, 0, 0 } };

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:b:", long_options, &option_index);

		if (key == -1) {
			break;
		}

		switch (key) {
			case 'a':
				alignmentsFilename = strdup(optarg);
				break;
            case 'b':
                sequencesFilename = strdup(optarg);
                break;
			default:
				return 1;
		}
	}

    FILE *sequencesFile = fopen(sequencesFilename, "r");
    struct List *seqs = constructEmptyList(0, NULL);
    struct List *seqLengths = constructEmptyList(0, free);
    struct List *headers = constructEmptyList(0, free);
    fastaRead(sequencesFile, seqs, seqLengths, headers);
    fclose(sequencesFile);

    FILE *alignmentsFile = fopen(alignmentsFilename, "r");
    int N = seqs->length;

    stHash *headerToID = stHash_construct();
    int64_t **sharedPositions = calloc(N, sizeof(int64_t*));
    for (int64_t i = 0; i < N; i++) {
        stHash_insert(headerToID, (void*) stHash_stringKey(headers->list[i]), (void*)i);
        sharedPositions[i] = (int64_t*) calloc(N, sizeof(int64_t));
    }
    struct PairwiseAlignment *alignment = NULL;
    while((alignment = cigarRead(alignmentsFile)) != NULL) {
        int64_t a = (int64_t) stHash_search(headerToID, (void*)stHash_stringKey(alignment->contig1));
        int64_t b = (int64_t) stHash_search(headerToID, (void*)stHash_stringKey(alignment->contig2));
        if (a == b) continue;
        int64_t seq1ID = (a < b) ? a : b;
        int64_t seq2ID = (a < b) ? b : a;
        
        int64_t matches = 0;
        for (int i = 0; i < alignment->operationList->length; i++) {
            struct AlignmentOperation *op = alignment->operationList->list[i];
            if (op->opType == PAIRWISE_MATCH) {
                matches += op->length;
            }
        }
        sharedPositions[seq1ID][seq2ID] += matches;
        free(alignment);
        
    }
    fclose(alignmentsFile);

    for (int64_t j = 0; j < N; j++) {
        for (int64_t i = 0; i < j; i++) {
            int64_t len_i = strlen(seqs->list[i]);
            int64_t len_j = strlen(seqs->list[j]);
            if (sharedPositions[i][j] == 0) continue;
            double p = 2*sharedPositions[i][j]/(double)(len_i + len_j);
            if (p <= 0.25) continue;
            double d =  (-3.0/4.0)*log(1.0 - (4.0/3.0)*(1.0 - p));
            fprintf(stdout, "%s %s %lf\n", (char*) headers->list[i], (char*) headers->list[j], d);
        }
    }
    
}