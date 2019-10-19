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
	stList *sequences = stList_construct();
	for (int i = 0; i < seqs->length; i++) {
		stList_append(sequences, seqs->list[i]);
	}

	stPinchThreadSet *graph = buildRepeatGraph(sequences, alignmentsFilename);

	if (gvizDebugFilename) {
		printBiedgedGraph(graph, gvizDebugFilename);
	}
	
	fprintf(stderr, "Graph has %ld blocks\n", stPinchThreadSet_getTotalBlockNumber(graph));
	stList *components = getOrdering(graph);
	for (int i = 0; i < stList_length(components); i++) {
		stList *component = stList_get(components, i);
		stList *path = heaviestPath(graph, component);
		stList *consensusSeq = traversePath(graph, path, sequences);
		fprintf(stdout, ">component_%d\n", i);
		printSeqList(consensusSeq);
		printf("\n");
	}

	stPinchThreadSet_destruct(graph);
}
