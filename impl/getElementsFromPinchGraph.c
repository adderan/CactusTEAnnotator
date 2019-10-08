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
	stList *components = getOrdering(threadSet);
	for (int i = 0; i < stList_length(components); i++) {
		stList *component = stList_get(components, i);
		stList *path = heaviestPath(graph, component);
		char *consensusSeq = traversePath(threadSet, path);
		fprintf(stdout, ">component_%d\n", i);
		fprintf(stdout, "%s\n", consensusSeq);
	}

	stPinchThreadSet_destruct(threadSet);
}
