#ifndef repeat_graph_h
#define repeat_graph_h
#include "stPinchGraphs.h"

#define _5PRIME 1
#define _3PRIME 0

#define SIDENAME(X) (X == 0) ? "3-prime" : "5-prime"

#define BLACK (void*)11
#define RED (void*)12
#define COLORNAME(X) (X == BLACK) ? "green" : "red"

#define OP(X) ((void*)X == BLACK ? RED : BLACK)

void pinchToGraphViz(stPinchThreadSet *threadSet, FILE *output);
stPinchThreadSet *buildRepeatGraph(stHash *sequences, char *alignmentsFilename);
stList *getOrdering(stPinchThreadSet *threadSet);
stList *getOrdering2(stPinchBlock *block);
stPinchBlock *getFirstBlock(stPinchThread *thread);
bool directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool startDirection);
void printBiedgedGraph(stPinchThreadSet *threadSet, char *gvizFilename);
stList *heaviestPath(stPinchThreadSet *graph, stList *ordering);
stList *traversePath(stPinchThreadSet *graph, stList *endsInPath, stHash *sequences);

#endif
