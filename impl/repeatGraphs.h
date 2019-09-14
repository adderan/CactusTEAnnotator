#ifndef repeat_graph_h
#define repeat_graph_h
#include "stPinchGraphs.h"

#define _5PRIME 1
#define _3PRIME 0

#define SIDENAME(X) (X == 0) ? "3-prime" : "5-prime"

void pinchToGraphViz(stPinchThreadSet *threadSet, FILE *output);
stPinchThreadSet *buildRepeatGraph(char *sequencesFilename, char *alignmentsFilename);
stPinchEnd *directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool direction);
stList *getOrdering(stPinchThreadSet *threadSet);
bool pinchCreatesCycle(stPinchSegment *seg1, stPinchSegment *seg2,
		bool orientation);


#endif
