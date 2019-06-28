#include <stdio.h>
#include "bioioC.h"
#include "pairwiseAlignment.h"

#include "stPinchIterator.h"
#include "stPinchGraphs.h"

int main(int argc, char **argv) {
	FILE *alignmentsFile = fopen(argv[1], "r");
	
	FILE *sequencesFile = fopen(argv[2], "r");
	struct List *seqs = constructEmptyList(0, NULL);
    struct List *seqLengths = constructEmptyList(0, free);
    struct List *headers = constructEmptyList(0, free);
	fastaRead(sequencesFile, seqs, seqLengths, headers);

	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	for (int i = 0; i < seqs->length; i++) {
		stPinchThreadSet_addThread(threadSet, stHash_stringKey(headers->list[i]), 0, strlen(seqs->list[i]));
	}

	//stPinchIterator *pinchIterator = stPinchIterator_constructFromFile(argv[1]);
	
	struct PairwiseAlignment *alignment = NULL;
	while((alignment = cigarRead(alignmentsFile)) != NULL) {

		stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, stHash_stringKey(alignment->contig1));
		stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, stHash_stringKey(alignment->contig2));
		int64_t length = alignment->end1 - alignment->start1;
		stPinchThread_pinch(thread1, thread2, alignment->start1, alignment->start2, length, alignment->strand1);

	}
	fprintf(stderr, "Total blocks: %ld\n", stPinchThreadSet_getTotalBlockNumber(threadSet));

	stPinchIterator_destruct(pinchIterator);
}
