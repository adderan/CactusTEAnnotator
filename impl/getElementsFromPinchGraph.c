#include "repeatGraphs.h"
#include "sonLib.h"
#include "bioioC.h"

int main(int argc, char **argv) {
	char *sequencesFilename = argv[1];
	char *alignmentsFilename = argv[2];
	char *gvizDebugFilename = argv[3];

	struct List *seqs = constructEmptyList(0, NULL);
	struct List *seqLengths = constructEmptyList(0, free);
	struct List *headers = constructEmptyList(0, free);
	FILE *sequencesFile = fopen(sequencesFilename, "r");
	fastaRead(sequencesFile, seqs, seqLengths, headers);
	fclose(sequencesFile);
	stHash *sequences = stHash_construct();
	int64_t threadName;
	for (int i = 0; i < seqs->length; i++) {
		sscanf(headers->list[i], "%ld", &threadName);
		stHash_insert(sequences, (void*) threadName, seqs->list[i]);
	}

	stPinchThreadSet *graph = buildRepeatGraph(sequences, alignmentsFilename);

	if (gvizDebugFilename) {
		//printBiedgedGraph(graph, gvizDebugFilename);
	}

	fprintf(stderr, "Graph has %ld blocks\n", stPinchThreadSet_getTotalBlockNumber(graph));

	assert(graphIsAcyclic(graph));

	stList *poGraph = getPartialOrderGraph(graph);
	assert(stPinchThreadSet_getTotalBlockNumber(graph) == stList_length(poGraph));
	
	int64_t i = 0;
	int64_t pathLength;
	do {
		stList *path = getHeaviestPath(poGraph);

		stList *consensusSequence = traversePath(graph, path, sequences);

		fprintf(stdout, ">consensus_%ld\n", i);
		pathLength = 0;
		for (int i = 0; i < stList_length(consensusSequence); i++) {
			char *subsequence = (char*)stList_get(consensusSequence, i);
			pathLength += strlen(subsequence);
			fprintf(stdout, "%s", subsequence);
		}
		fprintf(stdout, "\n");
		i++;
	}
	while (pathLength > 50);

	stPinchThreadSet_destruct(graph);
}
