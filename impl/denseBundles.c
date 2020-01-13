#include <stdio.h>
#include <stdlib.h>

#include "sonLib.h"
#include <getopt.h>

#include "lpo.h"
#include "seq_util.h"

#include "repeatGraphs.h"


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


/*
int *heaviestPath(LPOSequence_T *graph, int *pathLength, bool **seen) {

	int *containsPosition = calloc(sizeof(int), graph->nsource_seq);
	int *paths = calloc(sizeof(int), graph->length);
	int *score = calloc(sizeof(int), graph->length);
	LPOLetterSource_T *source;

	LPOLetter_T *seq = graph->letter;

	int bestScore = 0;
	int bestStart = -1;
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
		paths[i] = best_right;
		score[i] = right_score + right_overlap;

		if (score[i] > bestScore) {
			bestScore = score[i];
			bestStart = i;
		}
	}

	int *bestPath = calloc(sizeof(int), bestLength);
	int pos = bestStart;
	for (int i = 0; i < bestLength; i++) {
		bestPath[i] = pos;
		if (pos == bestEnd) break;
		pos = paths[pos];
	}
	*pathLength = bestLength;
	free(paths);
	free(score);
	free(containsPosition);

	return bestPath;

}
*/

void zeroPath(LPOSequence_T *graph, int *path, int pathLength) {
	LPOLetterSource_T *source;
	for (int i = 0; i < pathLength; i++) {
		source = &graph->letter[path[i]].source;
		do {
			graph->source_seq[source->iseq].weight = 0;
		}
		while ((source = source->more));
	}
}

//Get partial order graph out of POA data structure into
//my data structure
stList *lpoToPartialOrder(LPOSequence_T *graph) {
	LPOLetter_T *seq = graph->letter;

	stList *partialOrderGraph = stList_construct2(graph->length);
	for (int64_t i = 0; i < graph->length; i++) {
		PartialOrderNode *node = malloc(sizeof(PartialOrderNode));
		node->incomingNodes = stList_construct();
		node->incomingEdgeWeights = stList_construct();
		node->nodeID = i;
		node->length = 1;
		node->degree = 0;

		stList_set(partialOrderGraph, i, node);

		stHash *threads_i = stHash_construct();
		LPOLetterSource_T *source = &seq[i].source;
		while(source) {
			int64_t positionInSeq = (int64_t) source->ipos;
			int64_t seqID = (int64_t) source->iseq;
			stHash_insert(threads_i, (void*) seqID, (void*) positionInSeq);
			source = source->more;
			node->degree++;
		}
		assert(node->degree > 0);

		LPOLetterLink_T *left = NULL;
		for (left = &seq[i].left; (left) && (left->ipos >= 0); left= left->more) {
			PartialOrderNode *node_left = stList_get(partialOrderGraph, left->ipos);
			stList_append(node->incomingNodes, (void*) node_left->nodeID);

			source = &seq[left->ipos].source;
			int64_t edgeWeight = 0;
			while (source) {
				int64_t positionInSeq = (int64_t) source->ipos;
				int64_t seqID = (int64_t) source->iseq;
				if ((int64_t) stHash_search(threads_i, (void*) seqID) == positionInSeq + 1) {
					edgeWeight++;
				}
				source = source->more;
			}
			assert(edgeWeight > 0);
			stList_append(node->incomingEdgeWeights, (void*) edgeWeight);
		}

		stHash_destruct(threads_i);

	}
	return partialOrderGraph;
}

char *getConsensusSequence(stList *nodesInPath, LPOSequence_T *graph) {
	int64_t pathLength = stList_length(nodesInPath);
	char *path = calloc(pathLength + 1, sizeof(char));

	for (int64_t i = 0; i < stList_length(nodesInPath); i++) {
		PartialOrderNode *node = stList_get(nodesInPath, i);
		path[i] = graph->letter[node->nodeID].letter;
	}
	path[pathLength] = '\0';
	return path;
}

int main(int argc, char **argv) {
	char *lpoFilename = NULL;
	int64_t minConsensusLength = 50;
    while (1) {
        static struct option long_options[] = {
            { "lpo", required_argument, 0, 'a' }, 
			{ "minConsensusLength", required_argument, 0, 'b'},
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                lpoFilename = strdup(optarg);
                break;
			case 'b':
				sscanf(optarg, "%ld", &minConsensusLength);
				break;
            default:
                return 1;
        }
    }

	FILE *lpoFile = fopen(lpoFilename, "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

	int consensusNum = 0;

	stList *poGraph = lpoToPartialOrder(graph);
	char *consensusSeq = NULL;
	int64_t consensusLength = 0;
	do {
		stList *nodesInPath = getHeaviestPath(poGraph);
		printf(">consensus_%d\n", consensusNum);
		consensusSeq = getConsensusSequence(nodesInPath, graph);
		printf("%s\n", consensusSeq);
		consensusLength = strlen(consensusSeq);
		free(consensusSeq);
		consensusNum++;
	}
	while (consensusLength > minConsensusLength);
}
