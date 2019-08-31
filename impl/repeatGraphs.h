#ifndef repeat_graph_h
#define repeat_graph_h
#include "stPinchGraphs.h"

void pinchToGraphViz(stPinchThreadSet *threadSet, FILE *output);
stPinchThreadSet *buildRepeatGraph(char *sequencesFilename, char *alignmentsFilename);
bool directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool direction);
bool graphIsAcyclic(stPinchThreadSet *threadSet);

#endif
