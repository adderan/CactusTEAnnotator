#ifndef repeat_graph_h
#define repeat_graph_h
#include "stPinchGraphs.h"

#define _5PRIME 1
#define _3PRIME 0

#define SIDENAME(X) (X == 0) ? "3-prime" : "5-prime"

#define BLACK (void*)11;
#define RED (void*)12;
#define COLORNAME(X) (X == BLACK) ? "green" : "red"

#define OP(X) ((void*)X == BLACK ? RED : BLACK)

void pinchToGraphViz(stPinchThreadSet *threadSet, FILE *output);
stPinchThreadSet *buildRepeatGraph(char *sequencesFilename, char *alignmentsFilename);
bool directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool startDirection);
stList *getOrdering(stPinchThreadSet *threadSet, stHash *coloring);
void printBiedgedGraph(stPinchThreadSet *threadSet, char *gvizFilename);


#endif
