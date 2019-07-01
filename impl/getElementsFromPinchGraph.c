#include <stdio.h>
#include "bioioC.h"
#include "pairwiseAlignment.h"

#include "stPinchGraphs.h"
#include "stPinchIterator.h"

int main(int argc, char **argv) {
	
	FILE *sequencesFile = fopen(argv[2], "r");
	struct List *seqs = constructEmptyList(0, NULL);
    struct List *seqLengths = constructEmptyList(0, free);
    struct List *headers = constructEmptyList(0, free);
	fastaRead(sequencesFile, seqs, seqLengths, headers);

	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	int64_t threadName;
	for (int i = 0; i < seqs->length; i++) {
		sscanf(headers->list[i], "%ld", &threadName);
		stPinchThreadSet_addThread(threadSet, threadName, 0, strlen(seqs->list[i]));
	}
	
	stPinchIterator *pinchIterator = stPinchIterator_constructFromFile(argv[1]);
	stPinch *pinch = NULL;
	while((pinch = stPinchIterator_getNext(pinchIterator)) != NULL) {
		stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
		stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
		assert (thread1 != NULL);
		assert (thread2 != NULL);
		stPinchThread_pinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand);

	}
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


	fclose(sequencesFile);
	stPinchIterator_destruct(pinchIterator);
	stPinchThreadSet_destruct(threadSet);
}
