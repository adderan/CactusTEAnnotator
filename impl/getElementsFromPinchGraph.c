#include "repeatGraphs.h"
#include "sonLib.h"
#include "bioioC.h"


void printSeqList(stList *seqs) {
	for (int i = 0; i < stList_length(seqs); i++) {
		printf("%s", (char*) stList_get(seqs, i));
	}
}
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

	POGraph *poGraph = getPartialOrderGraph(graph);
	stList *path = heaviestPath(poGraph);

	fprintf(stderr, "Path length: %ld blocks\n", stList_length(path));


	stPinchThreadSet_destruct(graph);
}
