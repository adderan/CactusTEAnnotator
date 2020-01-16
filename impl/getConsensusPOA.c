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


char *getHeaviestPath(LPOSequence_T *graph, int64_t *nodeWeights, int64_t *bestPathScore) {

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

	*bestPathScore = 0;
	int64_t bestStart = 0;
	int64_t bestEnd = 0;
	int64_t bestPathLength = 0;
	for (int64_t start = 0; start < graph->length; start++) {
		stSet *seenThreads = stSet_construct();
		int64_t pathLength = 0;
		int64_t pathScore = 0;
		for (int64_t i = start; (i < graph->length) && (i != -1); i = direction[i]) {
			source = &seq[i].source;
			while(source) {
				char *threadName = graph->source_seq[source->iseq].name;
				if (!stSet_search(seenThreads, (void*) stHash_stringKey(threadName))) {
					stSet_insert(seenThreads, (void*) stHash_stringKey(threadName));
				}
				source = source->more;
			}
			pathLength++;
			int64_t nodeScore = 2*nodeWeights[i] - stSet_size(seenThreads);
			pathScore += nodeScore;

			if (pathScore > *bestPathScore) {
				bestStart = start;
				bestEnd = i;
				*bestPathScore = pathScore;
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
		nodeWeights[pos] = 0;
		pos = direction[pos];
		
	}

	free(direction);
	free(score);
	free(containsPosition);

	return path;
}

int main(int argc, char **argv) {
	char *lpoFilename = NULL;
	int64_t minConsensusScore = 100;
	char *namePrefix;
    while (1) {
        static struct option long_options[] = {
            { "lpo", required_argument, 0, 'a' }, 
			{ "minConsensusScore", required_argument, 0, 'b'},
			{ "namePrefix", required_argument, 0, 'c'},
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
				sscanf(optarg, "%ld", &minConsensusScore);
				break;
			case 'c':
				namePrefix = strdup(optarg);
				break;
            default:
                return 1;
        }
    }

	FILE *lpoFile = fopen(lpoFilename, "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

	int64_t *nodeWeights = calloc(graph->length, sizeof(int64_t));
	for (int64_t i = 0; i < graph->length; i++) {
		LPOLetterSource_T *source = &graph->letter[i].source;
		while (source) {
			nodeWeights[i]++;
			source = source->more;
		}
	}

	int consensusNum = 0;
	char *consensusSeq = NULL;
	int64_t consensusScore;
	while (true) {
		consensusSeq = getHeaviestPath(graph, nodeWeights, &consensusScore);
		if (consensusScore >= minConsensusScore) {
			printf(">%s_consensus_%d\n", namePrefix, consensusNum);
			printf("%s\n", consensusSeq);
		}
		free(consensusSeq);
		if (consensusScore < minConsensusScore) break;
		consensusNum++;
	}
}
