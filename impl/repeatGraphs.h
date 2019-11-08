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

typedef struct PONode {
	int64_t nodeID;
	int64_t *incomingNodes;
	int64_t nIncomingNodes;
	int64_t weight;
	bool orientation;
	void *data;
} PONode;

stPinchThreadSet *buildRepeatGraph(stHash *sequences, char *alignmentsFilename);
stList *getPartialOrderGraph(stPinchThreadSet *graph);
bool graphIsAcyclic(stPinchThreadSet *graph);
bool directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool startDirection);
//void printBiedgedGraph(stPinchThreadSet *threadSet, char *gvizFilename);
stList *heaviestPath(stList *poGraph);
stList *traversePath(stPinchThreadSet *graph, stList *endsInPath, stHash *sequences);

#endif
