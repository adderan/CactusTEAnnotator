#ifndef repeat_graph_h
#define repeat_graph_h
#include "stPinchGraphs.h"

#define MAX_BLOCK_DEGREE 100

#define _5PRIME 1
#define _3PRIME 0

#define SIDENAME(X) (X == 0) ? "3-prime" : "5-prime"

#define BLACK (void*)11
#define RED (void*)12
#define COLORNAME(X) (X == BLACK) ? "green" : "red"

#define OP(X) ((void*)X == BLACK ? RED : BLACK)


bool singleCopyFilterFn(stPinchSegment *seg1, stPinchSegment *seg2);
bool acyclicFilterFn(stPinchSegment *seg1, stPinchSegment *seg2);
bool graphIsAcyclic(stPinchThreadSet *graph);
bool directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool startDirection);
void printBiedgedGraph(stPinchThreadSet *threadSet, char *gvizFilename);
stList *getBlockOrdering(stPinchThreadSet *graph);
char *getConsensusSequence(stList *path, stHash *sequences);
stSortedSet *getConnectingThreads(stPinchEnd *end1, stPinchEnd *end2);
stList *tracebackHeaviestPath(stList *blockOrdering, int64_t *scores, int64_t *directions, int64_t *pathScore);
void getHeaviestPathScores(stList *blockOrdering, int64_t gapPenalty, int64_t *scores, int64_t *directions);
stPinchBlock *getHighestWeightBlock(stPinchThreadSet *graph);
stSortedSet *getThreads(stPinchSegment *segment);
#endif
