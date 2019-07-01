#include <stdio.h>
#include <stdlib.h>

#include "sonLib.h"

#include "lpo.h"
#include "seq_util.h"
#include <getopt.h>


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

int *dense_bundle(LPOSequence_T *graph, int *pathLength, double c) {

	int *containsPosition = calloc(sizeof(int), graph->nsource_seq);
	int *paths = calloc(sizeof(int), graph->length);
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
		paths[i] = best_right;
		score[i] = right_score + right_overlap;

		//printf("score = %d\n", score[i]);

	}

	double bestScore = 0.0;
	int bestStart, bestEnd, bestLength = 0;
	for (int start = 0; start < graph->length; start++) {
		int len = 0;
		for (int end = start; end >= 0; end = paths[end]) {
			len++;
			int weight = score[start] - score[end];
			double density = weight - c*len;
			if (density > bestScore) {
				bestScore = density;
				bestStart = start;
				bestEnd = end;
				bestLength = len;
			}
		}
	}
	if (bestScore <= 0) {
		return NULL;
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
	char *lpoFilename = NULL;
	bool iterate = false;
	double c = 1.1;
    while (1) {
        static struct option long_options[] = {
            { "lpo", required_argument, 0, 'a' }, 
			{ "iterate", no_argument, 0, 'b'},
			{ "lengthPenalty", required_argument, 0, 'c'},
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
				iterate = true;
                break;

			case 'c':
				sscanf(optarg, "%lf\n", &c);
				break;
            default:
                return 1;
        }
    }

	FILE *lpoFile = fopen(lpoFilename, "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

	int pathLength = 0;
	int *bestPath;
	int consensusNum = 0;
	if (iterate) {
		while (true) {
			//Keep extracting paths and zeroing out the weights
			//of all sequences in the best path until none are left
			bestPath = dense_bundle(graph, &pathLength, c);
			if (!bestPath) break;
			printf(">consensus_%d\n", consensusNum);
			for (int i = 0; i < pathLength; i++) {
				printf("%c", graph->letter[bestPath[i]].letter);
			}
			printf("\n");
			consensusNum++;
			zeroPath(graph, bestPath, pathLength);
		}
	}
	else {
		bestPath = dense_bundle(graph, &pathLength, c);
		for (int i = 0; i < pathLength; i++) {
			printf("%c", graph->letter[bestPath[i]].letter);
		}
		printf("\n");
	}
	
	free(bestPath);
}
