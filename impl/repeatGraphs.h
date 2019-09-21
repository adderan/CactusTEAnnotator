#ifndef repeat_graph_h
#define repeat_graph_h
#include "stPinchGraphs.h"

#define _5PRIME 1
#define _3PRIME 0

#define SIDENAME(X) (X == 0) ? "3-prime" : "5-prime"

#define BLACK (void*)1
#define RED (void*)2
#define COLORNAME(X) (X == (void*)1) ? "green" : "red"

#define OP(X) ((void*)X == BLACK ? RED : BLACK)

void pinchToGraphViz(stPinchThreadSet *threadSet, FILE *output);
stPinchThreadSet *buildRepeatGraph(char *sequencesFilename, char *alignmentsFilename);
bool directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool startDirection);
stList *getOrdering(stPinchThreadSet *threadSet, stHash *coloring);
bool pinchCreatesCycle(stPinchSegment *seg1, stPinchSegment *seg2,
		bool orientation);
void printBiedgedGraph(stPinchThreadSet *threadSet, 
		stHash *coloring, FILE *gvizFile);


#endif
