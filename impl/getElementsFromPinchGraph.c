#include <stdio.h>
#include "bioioC.h"
#include "pairwiseAlignment.h"

#include "stPinchGraphs.h"
#include "stPinchIterator.h"

stSortedSet *getThreads(stPinchBlock *block) {
	stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
	stSortedSet *threads = stSortedSet_construct();
	stPinchSegment *seg = NULL;
	while((seg = stPinchBlockIt_getNext(&blockIt)) != NULL) {
		stSortedSet_insert(threads, (void*) stPinchSegment_getName(seg));
	}
	return threads;
}

bool singleCopyFilterFn(stPinchSegment *seg1, stPinchSegment *seg2) {
	bool filter = false;
	stPinchBlock *block1 = stPinchSegment_getBlock(seg1);
	stPinchBlock *block2 = stPinchSegment_getBlock(seg2);

	if (!(block1 && block2)) return false;
	stSortedSet *threads1 = getThreads(block1);
	stSortedSet *threads2 = getThreads(block2);

	stSortedSet *intersect = stSortedSet_getIntersection(threads1, threads2);
	if (stSortedSet_size(intersect) > 0) filter = true;
	stSortedSet_destruct(threads1);
	stSortedSet_destruct(threads2);
	stSortedSet_destruct(intersect);
	return filter;
}

bool reachableOnThread(stPinchSegment *seg1, stPinchSegment *seg2, bool direction) {
	if (stPinchSegment_getThread(seg1) != stPinchSegment_getThread(seg2)) return false;
	bool beforeTarget = (stPinchSegment_getStart(seg1) <= (stPinchSegment_getStart(seg2) + stPinchSegment_getLength(seg2)));

	bool afterTarget = (stPinchSegment_getStart(seg2) <= stPinchSegment_getStart(seg1) + stPinchSegment_getLength(seg1)); 
	if (direction) {
		return beforeTarget;
	}
	return afterTarget;
}

/*

bool graphIsAcyclic(stPinchThreadSet *graph) {
	stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents(graph);
	stListIterator *componentsIt = stList_getIterator(adjacencyComponents);

	stList *component;
	while (component = stList_getNext(componentsIt)) {
		stSet *red = stSet_construct();
		stSet *black = stSet_construct();
		stPinchEnd *end = stList_peek(component);

		stList *stack = stList_construct();
		stPinchBlock *block;
		while(end) {
			block = stPinchEnd_getBlock(end);
			if (stSet_search(red, block)
		}
	}

	stList_destructIterator(componentsIt);
}

*/

stPinchEnd *getAdjacentEnd(stPinchSegment *segment, bool direction) {
	while (segment && (stPinchSegment_getBlock(segment) == NULL)) {
		segment = direction ? stPinchSegment_get3Prime(segment) : stPinchSegment_get5Prime(segment);
	}
	if (!segment) return NULL;
	stPinchBlock *block = stPinchSegment_getBlock(segment);
	if (!block) return NULL;

	bool orientation = !(direction ^ stPinchSegment_getBlockOrientation(segment));
	return stPinchEnd_construct(block, orientation);
}

//check if there is a directed walk from segment 1 to segment 2
bool directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool direction) {
	stPinchEnd *startEnd = getAdjacentEnd(seg1, direction);
	if (!startEnd) return false;
	printf("Starting at %d end of block at coordinate %ld.\n", stPinchEnd_getOrientation(startEnd), stPinchSegment_getStart(stPinchBlock_getFirst(stPinchEnd_getBlock(startEnd))));

	//Arriving at either of these ends means seg2 can
	//be traversed
	//First end encountered walking left from seg2
	stPinchEnd *leftTargetEnd = getAdjacentEnd(seg2, 1);
	//First end encountered walking right from seg2
	stPinchEnd *rightTargetEnd = getAdjacentEnd(seg2, 0);

	stPinchBlock *block = NULL;
	stSet *seen = stSet_construct();

	stList *stack = stList_construct();
	stPinchEnd *end = startEnd;
	stList_append(stack, end);
	bool reachable = false;
	while(stList_length(stack) > 0) {
		end = stList_pop(stack);
		block = stPinchEnd_getBlock(end);
		printf("Arrived at %d end of block at coordinate %ld.\n", stPinchEnd_getOrientation(end), stPinchSegment_getStart(stPinchBlock_getFirst(stPinchEnd_getBlock(end))));

		stPinchEnd otherEnd = stPinchEnd_constructStatic(block, !stPinchEnd_getOrientation(end));

		if ((leftTargetEnd && stPinchEnd_equalsFn(&otherEnd, leftTargetEnd)) || (rightTargetEnd && stPinchEnd_equalsFn(&otherEnd, rightTargetEnd))) {
			reachable = true;
			break;
		}
		if (!stSet_search(seen, block)) {
			//get the ends connected to the opposite side of this block
			//and add them to the stack
			stSet *adjEnds = stPinchEnd_getConnectedPinchEnds(&otherEnd);
			printf("Found %ld connected ends on %d end of block\n", stPinchEnd_getNumberOfConnectedPinchEnds(&otherEnd), stPinchEnd_getOrientation(&otherEnd));
			stSetIterator *adjEndIt = stSet_getIterator(adjEnds);
			stPinchEnd *adjEnd;
			while((adjEnd = stSet_getNext(adjEndIt)) != NULL) {
				stList_append(stack, adjEnd);
			}
			stSet_destructIterator(adjEndIt);
		}
	}

	stList_destruct(stack);
	return reachable;
}

void testDirectedWalk() {
	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	stPinchThreadSet_addThread(threadSet, 1, 0, 100);
	stPinchThreadSet_addThread(threadSet, 2, 0, 100);


	stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, 1);
	stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, 2);

	stPinchThread_pinch(thread1, thread2, 10, 10, 10, 1);

	stPinchSegment *seg1 = stPinchThread_getSegment(thread1, 25);
	stPinchSegment *seg2 = stPinchThread_getSegment(thread2, 25);

	//moving toward 5' end
	assert(!directedWalk(seg1, seg2, 0));
	//moving toward 3' end
	assert(!directedWalk(seg1, seg2, 1));

	stPinchThread_pinch(thread1, thread2, 50, 50, 10, 1);
	seg1 = stPinchThread_getSegment(thread1, 25);
	seg2 = stPinchThread_getSegment(thread2, 25);
	stPinchBlock *blockA = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 10));
	stPinchEnd blockALeftEnd = stPinchEnd_constructStatic(blockA, 1);
	assert(stPinchEnd_getNumberOfConnectedPinchEnds(&blockALeftEnd) == 0);
	assert(!directedWalk(seg1, seg2, 0));
	assert(!directedWalk(seg1, seg2, 1));


	stPinchThreadSet_destruct(threadSet);
	threadSet = stPinchThreadSet_construct();
	stPinchThreadSet_addThread(threadSet, 1, 0, 100);
	stPinchThreadSet_addThread(threadSet, 2, 0, 100);

	thread1 = stPinchThreadSet_getThread(threadSet, 1);
	thread2 = stPinchThreadSet_getThread(threadSet, 2);

	stPinchThread_pinch(thread1, thread2, 10, 10, 10, 0);

	seg1 = stPinchThread_getSegment(thread1, 25);
	seg2 = stPinchThread_getSegment(thread2, 25);

	//should be able to go backwards on thread 1, traverse
	//the reverse block, and then forwards on thread 2
	assert(directedWalk(seg1, seg2, 0));
	assert(!directedWalk(seg1, seg2, 1));

}

bool acyclicFilterFn(stPinchSegment *seg1, stPinchSegment *seg2) {
	if (singleCopyFilterFn(seg1, seg2)) return true;

	//dfs to check if block2 is reachable from block1
	if(directedWalk(seg1, seg2, 0) || directedWalk(seg1, seg2, 1)) return true;
	else return false;
}

void pinchToGraphViz(stPinchThreadSet *threadSet, FILE *output) {
	stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
	fprintf(output, "digraph {\n");

	stPinchBlock *block;
	while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
		stPinchEnd end = stPinchEnd_constructStatic(block, 0);
		stSet *adjacentEnds = stPinchEnd_getConnectedPinchEnds(&end);
		stSetIterator *endIterator = stSet_getIterator(adjacentEnds);
		stPinchEnd *adjacentEnd;
		while((adjacentEnd = stSet_getNext(endIterator)) != NULL) {

			fprintf(output, "\t%lu -> %lu\n", stHash_pointer(block), stHash_pointer(stPinchEnd_getBlock(adjacentEnd)));

		}

	}
	fprintf(output, "}\n");
}

bool isLeftStub(stPinchBlock *block) {
	stPinchEnd _3PrimeEnd = stPinchEnd_constructStatic(block, 0);
	return (stPinchEnd_getNumberOfConnectedPinchEnds(&_3PrimeEnd) == 0);
}

bool isRightStub(stPinchBlock *block) {
	stPinchEnd _5PrimeEnd = stPinchEnd_constructStatic(block, 1);
	return (stPinchEnd_getNumberOfConnectedPinchEnds(&_5PrimeEnd) == 0);
}

void highestWeightPath(stPinchThreadSet *threadSet) {
	stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);


	stPinchBlock *block = NULL;
	while((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
		if (isRightStub(block)) {
			printf("Found left stub %ld\n", stHash_pointer(block));
		}

	}
}

stPinchThreadSet *buildGraph(char *sequencesFilename, char *alignmentsFilename) {
	FILE *sequencesFile = fopen(sequencesFilename, "r");
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

	stPinchIterator *pinchIterator = stPinchIterator_constructFromFile(alignmentsFilename);
	stPinch *pinch = NULL;

	while((pinch = stPinchIterator_getNext(pinchIterator)) != NULL) {
		stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
		stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
		stPinchThread_filterPinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand, acyclicFilterFn);


	}
	fclose(sequencesFile);
	stPinchIterator_destruct(pinchIterator);
	return threadSet;
}

int main(int argc, char **argv) {
	/*
	stPinchThreadSet *threadSet = buildGraph(argv[1], argv[2]);

	if (argv[3] != NULL) {
		FILE *graphFile = fopen(argv[3], "w");
		pinchToGraphViz(threadSet, graphFile);
		fclose(graphFile);
	}


	stPinchThreadSet_destruct(threadSet);
	*/

	testDirectedWalk();
}
