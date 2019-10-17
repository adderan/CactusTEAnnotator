#include "repeatGraphs.h"

int main(int argc, char **argv) {
	char *sequencesFilename = argv[1];
	char *alignmentsFilename = argv[2];
	char *gvizDebugFilename = argv[3];

	stPinchThreadSet *graph = buildRepeatGraph(sequencesFilename, alignmentsFilename);

	if (gvizDebugFilename) {
		printBiedgedGraph(graph, gvizDebugFilename);
	}
	
	fprintf(stderr, "Graph has %ld blocks\n", stPinchThreadSet_getTotalBlockNumber(graph));
	stList *components = getOrdering(graph);
	for (int i = 0; i < stList_length(components); i++) {
		stList *component = stList_get(components, i);
		stList *path = heaviestPath(graph, component);
		//char *consensusSeq = traversePath(threadSet, path);
		fprintf(stdout, ">component_%d\n", i);
		fprintf(stdout, "%ld\n", stList_length(path));
	}

	stPinchThreadSet_destruct(graph);
}
