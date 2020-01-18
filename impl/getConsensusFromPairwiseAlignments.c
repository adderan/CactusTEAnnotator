#include "repeatGraphs.h"
#include "sonLib.h"
#include "bioioC.h"

#include "stPinchIterator.h"

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


	stPinchThreadSet *graph = stPinchThreadSet_construct();

	stHash *sequences = stHash_construct();
	
	for (int64_t i = 0; i < seqs->length; i++) {
		char *threadName = headers->list[i];
		int64_t threadID;
		sscanf(threadName, "%ld", &threadID);
		int64_t seqLength = strlen(seqs->list[i]);
		stHash_insert(sequences, (void*) threadID, seqs->list[i]);
		stPinchThreadSet_addThread(graph, threadID, 0, seqLength);
	}
	stPinchIterator *pinchIterator = 
		stPinchIterator_constructFromFile(alignmentsFilename);

	stPinch *pinch = NULL;
	while ((pinch = stPinchIterator_getNext(pinchIterator)) != NULL) {
		stPinchThread *thread1 = stPinchThreadSet_getThread(graph, pinch->name1);
		stPinchThread *thread2 = stPinchThreadSet_getThread(graph, pinch->name2);
		if (pinch->length < 10) continue;
		fprintf(stderr, "Pinching %ld to %ld\n", stPinchThread_getName(thread1), stPinchThread_getName(thread2));
		if (pinch->strand == 0) continue;
		stPinchThread_filterPinch(thread1, thread2, pinch->start1, pinch->start2, 
			pinch->length, pinch->strand, singleCopyFilterFn);
	}


	if (gvizDebugFilename) {
		printBiedgedGraph(graph, gvizDebugFilename);
	}

	fprintf(stderr, "Graph has %ld blocks\n", stPinchThreadSet_getTotalBlockNumber(graph));
	
	stList *blockOrdering = getBlockOrdering(graph);
	fprintf(stderr, "Ordering contains %ld blocks\n", stList_length(blockOrdering));

	stList *path = getHeaviestPath(blockOrdering, 1);
	fprintf(stderr, "Best path length %ld\n", stList_length(path));
	stPinchThreadSet_destruct(graph);
}
