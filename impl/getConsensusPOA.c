#include <stdio.h>
#include <stdlib.h>

#include "sonLib.h"
#include "bioioC.h"
#include <getopt.h>

#include "lpo.h"
#include "seq_util.h"


/*                  A
 *   A               \
 *    \               C-----G----A---C
 *     C ----T----G---C-----G----A---C--T----C---G
 *   G-C-----T        C-----G----A---C--
 *            \     /                   \
 *             \   /                     G---T---A  
 *               A
 *
 *  dense_bundle(G) = "CGAC"
 *
 *  heaviest_bundle(G) = "GCTGCGACTCG"
 *  
 *  heaviest_bundle is implemented by POA. It requires end-to-end
 *  traversal of the graph and does not care about the length of the
 *  consensus sequence produced, only optimizing the sum of weights
 *  of the path.
 *
 *  Dense bundle maximizes sum(score[x_i]) - c|x|
 *
 */


char *getHeaviestPath(LPOSequence_T *graph, stHash *seqLengths, double lengthPenalty) {

	int *containsPosition = calloc(sizeof(int), graph->nsource_seq);
	int *direction = calloc(sizeof(int), graph->length);
	memset(direction, -1, graph->length);
	int *score = calloc(sizeof(int), graph->length);
	LPOLetterSource_T *source;

	LPOLetter_T *seq = graph->letter;

	for (int i = graph->length - 1; i >= 0; i--) {
		source = &seq[i].source;
		memset(containsPosition, 0, graph->nsource_seq * sizeof(int));
		do {
			if (graph->source_seq[source->iseq].weight > 0) {
				containsPosition[source->iseq] = source->ipos;
			}
		}
		while ((source = source->more));

		int right_score = 0;
		int right_overlap = 0;
		int best_right = -1;

		LPOLetterLink_T *right = NULL;
		for (right = &seq[i].right; (right) && (right->ipos >= 0); right= right->more) {
			int overlap = 0;
			source = &seq[right->ipos].source;
			do {
				//check if nodes seq[i] and this right neighbor
				//are directly adjacent in this sequence
				if (containsPosition[source->iseq] + 1 == source->ipos) {
					overlap += graph->source_seq[source->iseq].weight;
				}
			}
			while ((source = source->more));

			if ((overlap > right_overlap) ||
					(overlap == right_overlap && score[right->ipos] > right_score)) {
				right_overlap = overlap;
				right_score = score[right->ipos];
				best_right = right->ipos;
			}
		}
		direction[i] = best_right;
		score[i] = right_score + right_overlap;

	}

	int64_t bestStart = 0;
	int64_t bestEnd = 0;
	int64_t bestPathLength = 0;
	double bestScore = 0.0;
	for (int64_t start = 0; start < graph->length; start++) {
		stSet *seenThreads = stSet_construct();
		int64_t totalPositions = 0;
		int64_t pathLength = 0;
		int64_t pathWeight = 0;
		for (int64_t i = start; (i < graph->length) && (i != -1); i = direction[i]) {
			source = &seq[i].source;
			int64_t nodeWeight = 0;
			while(source) {
				nodeWeight += graph->source_seq[source->iseq].weight;
				if (!stSet_search(seenThreads, (void*) (int64_t) source->iseq)) {
					stSet_insert(seenThreads, (void*) (int64_t) source->iseq);
					char *threadName = graph->source_seq[source->iseq].name;

					int64_t threadLength = (int64_t) stHash_search(seqLengths, (void*)stHash_stringKey(threadName));
					totalPositions += threadLength;
				}
				source = source->more;
			}

			pathLength++;
			pathWeight += nodeWeight;
			

			double pathScore = (double)pathWeight/(double)totalPositions;
			if (pathScore > bestScore) {
				bestStart = start;
				bestEnd = i;
				bestScore = pathScore;
				bestPathLength = pathLength;
			}


		}
		stSet_destruct(seenThreads);
	}


	char *path = calloc(bestPathLength + 1, sizeof(char));

	int64_t i = 0;
	int64_t pos = bestStart;
	while ((pos <= bestEnd) && pos != -1) {
		path[i] = seq[pos].letter;
		i++;

		source = &seq[pos].source;
		while (source != NULL) {
			graph->source_seq[source->iseq].weight = 0;
			source = source->more;
		}

		pos = direction[pos];
		
	}
	path[bestPathLength] = '\0';

	free(direction);
	free(score);
	free(containsPosition);

	return path;
}

int main(int argc, char **argv) {
	char *lpoFilename = NULL;
	char *sequencesFilename = NULL;
	int64_t minConsensusLength = 5;
    while (1) {
        static struct option long_options[] = {
            { "lpo", required_argument, 0, 'a' }, 
			{ "sequences", required_argument, 0, 'b'},
			{ "minConsensusLength", required_argument, 0, 'c'},
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                lpoFilename = strdup(optarg);
                break;
			case 'b':
				sequencesFilename = strdup(optarg);
				break;
			case 'c':
				sscanf(optarg, "%ld", &minConsensusLength);
				break;
            default:
                return 1;
        }
    }

	FILE *lpoFile = fopen(lpoFilename, "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

	FILE *sequencesFile = fopen(sequencesFilename, "r");
    struct List *seqsList = constructEmptyList(0, NULL);
    struct List *seqLengthsList = constructEmptyList(0, free);
    struct List *headersList = constructEmptyList(0, free);
    fastaRead(sequencesFile, seqsList, seqLengthsList, headersList);
    fclose(sequencesFile);

	stHash *seqLengths = stHash_construct();
	for (int64_t i = 0; i < seqsList->length; i++) {
		stHash_insert(seqLengths, (void*)stHash_stringKey(headersList->list[i]), (void*) strlen(seqsList->list[i]));
	}

	int consensusNum = 0;

	char *consensusSeq = NULL;
	while (true) {
		consensusSeq = getHeaviestPath(graph, seqLengths, 1.0);
		if (strlen(consensusSeq) == 0) break;

		printf(">consensus_%d\n", consensusNum);
		printf("%s\n", consensusSeq);
		free(consensusSeq);
		consensusNum++;

	}
}
