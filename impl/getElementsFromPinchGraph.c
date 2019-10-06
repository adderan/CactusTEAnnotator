#include "repeatGraphs.h"

int main(int argc, char **argv) {
	char *sequencesFilename = argv[1];
	char *alignmentsFilename = argv[2];
	char *gvizDebugFilename = argv[3];

	stPinchThreadSet *threadSet = buildRepeatGraph(sequencesFilename, alignmentsFilename);

	if (gvizDebugFilename) {
		printBiedgedGraph(threadSet, gvizDebugFilename);
	}
	
	fprintf(stderr, "Graph has %ld blocks\n", stPinchThreadSet_getTotalBlockNumber(threadSet));
	stList *orderings = getOrdering(threadSet);
	stListIterator *it = stList_getIterator(orderings);
	for (int i = 0; i < stList_length(components); i++) {
		stList *component = stList_search(components, i);
		char *path = getConsensusPath(threadSet, componentOrdering);
		fprintf(stdout, ">component_%d\n", i);
		fprintf(stdout, "%s\n", path);
	}

	stPinchThreadSet_destruct(threadSet);
}
