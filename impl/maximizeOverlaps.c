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


int main(int argc, char **argv) {

    FILE *lpoFile = fopen(argv[1], "r");
    LPOSequence_T *graph = read_lpo(lpoFile);
    fclose(lpoFile);

    LPOLetter_T *seq = graph->letter;

    int *path = malloc(sizeof(int)*graph->length);


    int **score = malloc(graph->length, sizeof(int*));
    for (int i = 0; i < graph->length; i++) {
        score[i] = calloc(graph->length, sizeof(int));
    }


    int bestScore = -999999;
    int bestNode = -1;

    fprintf(stderr, "Graph length = %d\n", graph->length);


    int *seenThreads = (int*) calloc(sizeof(int), graph->nsource_seq);


    for (int i = graph->length - 1; i >= 0; i--) {

        int newThreads = 0;
        for (LPOLetterSource_T *source = &seq[i].source; source!= NULL; source = source->more) {
            if (seenThreads[source->iseq] == 0) {
                newThreads++;
                seenThreads[source->iseq] = 1;
            }
        }
        printf("Found %d new threads at node %d\n", newThreads, i);

        int maxWeight = 0;
        int bestEdge = -1;
        int rightScore = 0;
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
        score[i] = rightScore + maxWeight - newThreads;
        if (score[i] > bestScore) {
            bestScore = score[i];
            bestNode = i;
        }
        printf("Score of node %d: %d\n", i, score[i]);
    }

    //traceback the best path
    int *bestPath = (int*) calloc(sizeof(int), graph->length);

    //start from the highest scoring node in the graph
    int i = bestNode;
    int pathLength = 0;
    while(i >= 0) {
        bestPath[pathLength++] = i;
        i = path[i];
    }

    for (int i = 0; i < pathLength; i++) {
        printf("%c", seq[bestPath[i]].letter);
    }
    printf("\n");

}
