#include <stdio.h>
#include "bioioC.h"
#include "pairwiseAlignment.h"

#include "stPinchGraphs.h"
#include "stPinchIterator.h"

#include "repeatGraphs.h"


bool pinchIsNegative;

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

stList *getComponentOrdering2(stPinchEnd *startEnd) {
	stList *stack = stList_construct();
	stList_append(stack, startEnd);

	stSet *red = stSet_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, (void(*)(void *)) stPinchEnd_destruct);
	stSet *black = stSet_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, (void(*)(void *)) stPinchEnd_destruct);

	stList *ordering = stList_construct();

	stPinchBlock *block;
	stPinchEnd *end;
	while (stList_length(stack) > 0) {
		end = stList_pop(stack);
		fprintf(stderr, "Arrived at end %d of block starting at %ld\n", stPinchEnd_getOrientation(end), stPinchSegment_getStart(stPinchBlock_getFirst(stPinchEnd_getBlock(end))));
		block = stPinchEnd_getBlock(end);

		bool seenBlock = false;
		if (stSet_search(black, end)) {
			fprintf(stderr, "End is black\n");
			seenBlock = true;
		}
		else if (stSet_search(red, end)) {
			fprintf(stderr, "Encountered cycle\n");
			stList_destruct(ordering);
			stSet_destruct(red);
			stSet_destruct(black);
			return NULL;

		}
		else {
			fprintf(stderr, "End is unlabeled, labeling it black\n");
			stSet_insert(black, end);
			stList_append(ordering, end);
		}

		stPinchEnd *oppositeEnd = stPinchEnd_construct(block, !stPinchEnd_getOrientation(end));
		if (!stSet_search(red, oppositeEnd)) {
			stSet_insert(red, oppositeEnd);
		}
		if (!seenBlock) {
			stSet *adjacentEnds = stPinchEnd_getConnectedPinchEnds(oppositeEnd);
			fprintf(stderr, "Found %ld ends adjacent to %d side\n", stSet_size(adjacentEnds), stPinchEnd_getOrientation(oppositeEnd));
			stSetIterator *adjEndIt = stSet_getIterator(adjacentEnds);
			stPinchEnd *adjEnd;
			while((adjEnd = stSet_getNext(adjEndIt))) {
				stList_append(stack, adjEnd);
			}
			stSet_destructIterator(adjEndIt);
		}
	}
	stSet_destruct(red);
	stSet_destruct(black);
	return ordering;
}

stList *getComponentOrdering(stPinchBlock *block) {
	stPinchEnd *leftStart = stPinchEnd_construct(block, 0);
	stPinchEnd *rightStart = stPinchEnd_construct(block, 1);
	stList *leftOrdering = getComponentOrdering2(leftStart);
	stList *rightOrdering = getComponentOrdering2(rightStart);
	if (!(leftOrdering && rightOrdering)) return NULL;
	stList *ordering = stList_construct();
	stList_appendAll(ordering, leftOrdering);
	stList_appendAll(ordering, rightOrdering);
	return ordering;
}

/*Get an ordering of the nodes in the graph, such that each node
 * is preceded in the ordering by all nodes from which it is accessible
 * by a directed walk. Return NULL if the graph is not acyclic. */
stList *getOrdering(stPinchThreadSet *graph) {
	stSortedSet *threadComponents = stPinchThreadSet_getThreadComponents(graph);
	fprintf(stderr, "Found %ld components\n", stSortedSet_size(threadComponents));
	stList *ordering = stList_construct();

	stSortedSetIterator *componentsIt = stSortedSet_getIterator(threadComponents);
	stList *component;
	while ((component = stSortedSet_getNext(componentsIt))) {

		//get any block in the component
		stPinchSegment *segment = stPinchThread_getFirst(stList_peek(component));
		while (segment && (stPinchSegment_getBlock(segment) == NULL)) segment = stPinchSegment_get3Prime(segment);

		stPinchBlock *block = stPinchSegment_getBlock(segment);
		stList *componentOrdering = getComponentOrdering(block);
		if (!componentOrdering) return NULL;
		stList_appendAll(ordering, componentOrdering);

	}
	stSortedSet_destructIterator(componentsIt);
	return ordering;
}

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

bool directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool direction) {
	stPinchEnd *startEnd = getAdjacentEnd(seg1, direction);
	if (!startEnd) return false;
	fprintf(stderr, "Starting at %d end of block at coordinate %ld.\n", stPinchEnd_getOrientation(startEnd), stPinchSegment_getStart(stPinchBlock_getFirst(stPinchEnd_getBlock(startEnd))));

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
		fprintf(stderr, "Arrived at %d end of block at coordinate %ld.\n", stPinchEnd_getOrientation(end), stPinchSegment_getStart(stPinchBlock_getFirst(stPinchEnd_getBlock(end))));

		stPinchEnd otherEnd = stPinchEnd_constructStatic(block, !stPinchEnd_getOrientation(end));

		if ((leftTargetEnd && stPinchEnd_equalsFn(&otherEnd, leftTargetEnd)) || (rightTargetEnd && stPinchEnd_equalsFn(&otherEnd, rightTargetEnd))) {
			reachable = true;
			break;
		}
		if (!stSet_search(seen, block)) {
			stSet_insert(seen, block);
			//get the ends connected to the opposite side of this block
			//and add them to the stack
			stSet *adjEnds = stPinchEnd_getConnectedPinchEnds(&otherEnd);
			fprintf(stderr, "Found %ld connected ends on %d end of block\n", stPinchEnd_getNumberOfConnectedPinchEnds(&otherEnd), stPinchEnd_getOrientation(&otherEnd));
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

bool acyclicFilterFn(stPinchSegment *seg1, stPinchSegment *seg2) {
	if (singleCopyFilterFn(seg1, seg2)) return true;

	if (pinchIsNegative) {
		return directedWalk(seg1, seg2, 0);
	}
	return directedWalk(seg1, seg2, 1);
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

stPinchThreadSet *buildRepeatGraph(char *sequencesFilename, char *alignmentsFilename) {
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

		//need to store this in a global because filterPinch doesn't
		//provide the orientation to the filter function
		pinchIsNegative = pinch->strand;

		stPinchThread_filterPinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand, acyclicFilterFn);


	}
	fclose(sequencesFile);
	stPinchIterator_destruct(pinchIterator);
	return threadSet;
}


