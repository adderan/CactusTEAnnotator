#include <stdio.h>
#include <stdlib.h>

#include "sonLib.h"

#include "lpo.h"
#include "seq_util.h"


int edgeWeight(LPOSequence_T *graph, LPOLetter_T *left, LPOLetter_T *right) {
	int weight = 0;
	int *containsPosition = (int*) calloc(sizeof(int), graph->nsource_seq);

	LPOLetterSource_T *leftSource = &left->source;
	LPOLetterSource_T *rightSource = &right->source;
	do {
		containsPosition[leftSource->iseq] = leftSource->ipos+1;
	}
	while (leftSource = leftSource->more);

	do {
		if (containsPosition[rightSource->iseq] == rightSource->ipos) {
			weight++;
		}
	}
	while (rightSource = rightSource->more);

	free(containsPosition);
	return weight;
}


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


int main(int argc, char **argv) {

	FILE *lpoFile = fopen(argv[1], "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

	LPOLetter_T *seq = graph->letter;

	int *path = malloc(sizeof(int)*graph->length);


	int *score = (int*)calloc(graph->length, sizeof(int));

	int *nThreads = (int*)calloc(graph->length, sizeof(int));

	int **seenThreads = (int**)calloc(graph->length, sizeof(int*));
	for (int i = 0; i < graph->length; i++) {
		seenThreads[i] = (int*)calloc(graph->nsource_seq, sizeof(int));
	}

	int bestScore = -999999;
	int bestNode = -1;

	fprintf(stderr, "Graph length = %d\n", graph->length);

	for (int i = graph->length - 1; i >= 0; i--) {
		int maxWeight = 0;
		int bestEdge = -1;
		int rightScore = 0;
		for (int j = 0; j < graph->nsource_seq; j++) {
			if (i+1 < graph->length) {
				seenThreads[i][j] = seenThreads[i+1][j];
			}
		}
		for (LPOLetterSource_T *source = &seq[i].source; source != NULL; source = source->more) {
			seenThreads[i][source->iseq] = 1;
		}
		for (int j = 0; j < graph->nsource_seq; j++) {
			nThreads[i] += seenThreads[i][j];
		}
		fprintf(stderr, "Seen %d threads\n", nThreads[i]);
		for (LPOLetterLink_T *right = &seq[i].right; right && right->ipos>=0; 
				right=right->more) {

			int weight = edgeWeight(graph, &seq[i], &seq[right->ipos]);
			if ((weight > maxWeight) || ((weight == maxWeight) && (score[right->ipos] > rightScore))) {
				maxWeight = weight;
				bestEdge = right->ipos;
				rightScore = score[right->ipos];
			}

		}
		path[i] = bestEdge;
		score[i] = rightScore + maxWeight;
		if (score[i] > bestScore) {
			bestScore = score[i];
			bestNode = i;
		}
		printf("Score of node %d: %d\n", i, score[i]);
	}

	//traceback the best path
	int *bestPath = (int*) calloc(sizeof(int), graph->length);

	int bestStart = 0;
	int bestEnd = 0;
	double highestDensity = 0.0;
	for (int start = 0; start < graph->length; start++) {
		int currentNode = path[start];
		for (int end = start; end < graph->length; end = path[currentNode++]) {
			int newThreads = nThreads[end] - nThreads[start];
			if (newThreads == 0) {
				newThreads = 1;
			}
			double density = (double)(score[start] - score[end])/newThreads;
			//printf("score start = %d\n", score[start]);
			//printf("score end = %d\n", score[end]);
			//printf("seen threads start = %d\n", nThreads[start]);
			//printf("seen threads end = %d\n", nThreads[end]);
			//printf("Density = %f\n", density);
			if (density > highestDensity) {
				highestDensity = density;
				bestStart = start;
				bestEnd = end;
			}
		}
	}
	printf("Start = %d\n", bestStart);
	printf("End = %d\n", bestEnd);
	printf("Highest density = %f\n", highestDensity);

	int i = bestStart;
	int pathLength = 0;
	while(i >= 0 && i != bestEnd) {
		bestPath[pathLength++] = i;
		i = path[i];
	}

	for (int i = 0; i < pathLength; i++) {
		printf("%c", seq[bestPath[i]].letter);
	}
	printf("\n");

}
