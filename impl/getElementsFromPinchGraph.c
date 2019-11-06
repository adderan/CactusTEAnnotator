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
<<<<<<< HEAD
	fprintf(stderr, "Total blocks: %ld\n", stPinchThreadSet_getTotalBlockNumber(threadSet));

	fprintf(stderr, "Number of components: %ld\n", stSortedSet_size(stPinchThreadSet_getThreadComponents(threadSet)));

	stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
	//stPinchBlock *block;
	//while((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
	//	printf("Block degree: %ld\n", stPinchBlock_getDegree(block));
	//}
	
	stPinchBlock *firstBlock = stPinchThreadSetBlockIt_getNext(&blockIt);
	stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(firstBlock);
	stPinchSegment *firstSegment = stPinchBlockIt_getNext(&segmentIt);
	printf("First segment start %ld\n", stPinchSegment_getStart(firstSegment));
	stPinchEnd firstEnd = stPinchEnd_constructStatic(firstBlock, 1);
	stPinchEnd *end = &firstEnd;
	stPinchBlock *block;
	stSet *seen = stSet_construct();
	while(true) {
		block = stPinchEnd_getBlock(end);
		printf("Arrived at %d end of block of length %ld, degree %ld\n", stPinchEnd_getOrientation(end), stPinchBlock_getLength(block), stPinchBlock_getDegree(block));
		if (stSet_search(seen, block)) {
			printf("Traversed cycle");
			break;
		}
		stSet_insert(seen, block);

		stSet *adjacentEnds = stPinchEnd_getConnectedPinchEnds(end);
		printf("Found %ld connected blocks\n", stSet_size(adjacentEnds));
		stPinchEnd *nextEnd = stSet_peek(adjacentEnds);
		if (nextEnd == NULL) {
			printf("No adjacent ends\n");
			break;
		}
		printf("Number of subsequences: %ld\n", stList_length(stPinchEnd_getSubSequenceLengthsConnectingEnds(end, nextEnd)));
		end = nextEnd;

	}

=======
>>>>>>> pinch-graphs-redo-ordering

	fprintf(stderr, "Graph has %ld blocks\n", stPinchThreadSet_getTotalBlockNumber(graph));

	assert(graphIsAcyclic(graph));

	stList *poGraph = getPartialOrderGraph(graph);
	assert(stPinchThreadSet_getTotalBlockNumber(graph) == stList_length(poGraph));
	stList *path = heaviestPath(poGraph);

	fprintf(stderr, "Path length: %ld blocks\n", stList_length(path));


	stPinchThreadSet_destruct(graph);
}
