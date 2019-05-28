#include <stdio.h>
#include <stdlib.h>

#include "sonLib.h"

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
 *
 */

int *getDenseBundle(LPOSequence_T *graph, int *pathLength) {

	int *containsPosition = calloc(sizeof(int), graph->nsource_seq);
	int *paths = calloc(sizeof(int), graph->length);
	int *score = calloc(sizeof(int), graph->length);
	LPOLetterSource_T *source;

	LPOLetter_T *seq = graph->letter;
	int ibest = -1;
	int best_score = -1;

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
		if (score[i] > best_score) {
			ibest = i;
			best_score = score[i];
		}

	}
	if (best_score <= 0) {
		return NULL;
	}

	int *bestPath = calloc(sizeof(int), graph->length);
	while (ibest >= 0) {
		bestPath[*pathLength] = ibest;
		*pathLength = *pathLength + 1;
		ibest = paths[ibest];
	}
	free(paths);
	free(score);
	free(containsPosition);

	return bestPath;

}

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


int main(int argc, char **argv) {

	FILE *lpoFile = fopen(argv[1], "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

	int pathLength = 0;
	int *bestPath;
	int consensusNum = 0;
	while (true) {
		//Keep extracting paths and zeroing out the weights
		//of all sequences in the best path until none are left
		bestPath = getDenseBundle(graph, &pathLength);
		if (!bestPath) break;
		printf(">consensus_%d\n", consensusNum);
		for (int i = 0; i < pathLength; i++) {
			printf("%c", graph->letter[bestPath[i]].letter);
		}
		printf("\n");
		consensusNum++;
		zeroPath(graph, bestPath, pathLength);
	}
	
	free(bestPath);
}
