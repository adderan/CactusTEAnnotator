#include "repeatGraphs.h"

int main(int argc, char **argv) {
	stPinchThreadSet *threadSet = buildRepeatGraph(argv[1], argv[2]);

	if (argv[3] != NULL) {
		FILE *graphFile = fopen(argv[3], "w");
		pinchToGraphViz(threadSet, graphFile);
		fclose(graphFile);
	}
	/*
	stList *ordering = getOrdering(threadSet);
	assert(ordering);
	for (int i = 0; i < stList_length(ordering); i++) {
		fprintf(stderr, "%p ", stList_get(ordering, i));
	}
	fprintf(stderr, "\n");
	*/

	stPinchThreadSet_destruct(threadSet);
}
