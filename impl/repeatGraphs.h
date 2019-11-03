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

typdedef struct POANode {
	int64_t nodeID;
	int64_t incomingNodes[];
	int64_t weight;
	stPinchBlock *block;
} POANode;

void pinchToGraphViz(stPinchThreadSet *threadSet, FILE *output);
stPinchThreadSet *buildRepeatGraph(stHash *sequences, char *alignmentsFilename);
stList *getDAG(stPinchThreadSet *graph);
bool graphIsAcyclic(stPinchThreadSet *graph);
bool directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool startDirection);
void printBiedgedGraph(stPinchThreadSet *threadSet, char *gvizFilename);
stList *heaviestPath(stPinchThreadSet *graph, stList *ordering);
stList *traversePath(stPinchThreadSet *graph, stList *endsInPath, stHash *sequences);

#endif
