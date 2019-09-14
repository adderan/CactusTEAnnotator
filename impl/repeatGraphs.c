#include <stdio.h>
#include "bioioC.h"
#include "pairwiseAlignment.h"

#include "stPinchGraphs.h"
#include "stPinchIterator.h"

#include "repeatGraphs.h"

bool pinchOrientation;

char *endInfo(stPinchEnd *end) {
	stPinchBlock *block = stPinchEnd_getBlock(end);
	stPinchSegment *seg = stPinchBlock_getFirst(block);
	int segStart = stPinchSegment_getStart(seg);
	int segEnd = stPinchSegment_getStart(seg) + stPinchSegment_getLength(seg);
	char *side = SIDENAME(stPinchEnd_getOrientation(end));
	char *info = malloc(sizeof(char)*100);
	sprintf(info, "%s end of block from %d-%d", side, segStart, segEnd);
	return info;
}

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

/*If the graph is not acyclic, return NULL. If the graph is acyclic, 
 * return an ordering of the nodes such that all nodes reachable 
 * from any given node occur after it in the ordering.*/
#define BLACK (void*)1
#define RED (void*)2

#define OP(X) ((void*)X == BLACK ? RED : BLACK)
stList *getComponentOrdering(stPinchBlock *startBlock) {
	stList *stack = stList_construct();
	stPinchEnd *leftEnd = stPinchEnd_construct(startBlock, _5PRIME);
	stPinchEnd *rightEnd = stPinchEnd_construct(startBlock, _3PRIME);
	stList_append(stack, leftEnd);
	stList_append(stack, rightEnd);

	stSet *seen = stSet_construct();

	stHash *color = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, (void(*)(void*)) stPinchEnd_destruct, NULL);

	stHash_insert(color, leftEnd, BLACK);
	stHash_insert(color, rightEnd, RED);

	stList *ordering = stList_construct();

	stPinchBlock *block;
	stPinchEnd *end;
	while (stList_length(stack) > 0) {
		end = stList_pop(stack);
		fprintf(stderr, "Arrived at %s\n", endInfo(end));
		block = stPinchEnd_getBlock(end);

		stPinchEnd *oppositeEnd = stPinchEnd_construct(block, !stPinchEnd_getOrientation(end));

		assert(stHash_search(color, end));
		if (stHash_search(color, end) ==
				stHash_search(color, oppositeEnd)) {
			//encountered cycle
			fprintf(stderr, "Encountered cycle.\n");
			stList_destruct(stack);
			stHash_destruct(color);
			stList_destruct(ordering);
			return NULL;
		}
		else if (!stHash_search(color, oppositeEnd)) {
			stHash_insert(color, oppositeEnd, OP(stHash_search(color, end)));
		}
		if (!stSet_search(seen, block)) {
			stSet_insert(seen, block);
			stSet *nearAdjacentEnds = stPinchEnd_getConnectedPinchEnds(end);
			fprintf(stderr, "Found %ld connected ends on near side\n", stSet_size(nearAdjacentEnds));
			stSetIterator *adjEndIt = stSet_getIterator(nearAdjacentEnds);
			stPinchEnd *adjEnd;
			while((adjEnd = stSet_getNext(adjEndIt))) {
				stHash_insert(color, adjEnd, OP(stHash_search(color, end)));
				stList_append(stack, adjEnd);
			}

			stSet_destructIterator(adjEndIt);
			stSet *farAdjacentEnds = stPinchEnd_getConnectedPinchEnds(oppositeEnd);
			fprintf(stderr, "Found %ld connected ends on far side.\n", stSet_size(farAdjacentEnds));
			adjEndIt = stSet_getIterator(farAdjacentEnds);
			while((adjEnd = stSet_getNext(adjEndIt))) {
				stHash_insert(color, adjEnd, OP(stHash_search(color, oppositeEnd)));
				stList_append(stack, adjEnd);
			}
			stSet_destructIterator(adjEndIt);
		}
	}
	stHash_destruct(color);
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
		if (stList_length(component) == 1) {
			//only one thread with no blocks, so nothing to order
			continue;
		}

		//get any block in the component
		stPinchSegment *segment = stPinchThread_getFirst(stList_peek(component));
		while (segment && (stPinchSegment_getBlock(segment) == NULL)) segment = stPinchSegment_get3Prime(segment);
		if (segment == NULL) {
			//try other direction
			segment = stPinchThread_getFirst(stList_peek(component));
			while (segment && (stPinchSegment_getBlock(segment) == NULL)) {
				segment = stPinchSegment_get3Prime(segment);
			}
		}
		//The case of a single thread with no blocks has already
		//been handled, so there must be a segment contained
		//in a block
		assert(segment);

		stPinchBlock *block = stPinchSegment_getBlock(segment);
		stList *componentOrdering = getComponentOrdering(block);
		if (!componentOrdering) return NULL;
		stList_appendAll(ordering, componentOrdering);

	}
	stSortedSet_destructIterator(componentsIt);
	return ordering;
}

stPinchEnd *getAdjacentEnd(stPinchSegment *segment, bool direction) {
	if (stPinchSegment_getBlock(segment)) {
		return stPinchEnd_construct(stPinchSegment_getBlock(segment), direction);
	}
	while (segment && (stPinchSegment_getBlock(segment) == NULL)) {
		if (direction == _3PRIME) {
			segment = stPinchSegment_get3Prime(segment);
		}
		else {
			assert(direction == _5PRIME);
			segment = stPinchSegment_get5Prime(segment);
		}
	}
	if (!segment) return NULL;
	stPinchBlock *block = stPinchSegment_getBlock(segment);
	if (!block) return NULL;

	bool orientation = direction ^ stPinchSegment_getBlockOrientation(segment);
	return stPinchEnd_construct(block, orientation);
}

stPinchEnd *directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool direction) {
	stPinchEnd *startEnd = getAdjacentEnd(seg1, direction);
	if (!startEnd) return false;

	//Arriving at either of these ends means seg2 can
	//be traversed
	//First end encountered walking left from seg2
	stPinchEnd *leftTargetEnd = getAdjacentEnd(seg2, _5PRIME);
 
	//First end encountered walking right from seg2
	stPinchEnd *rightTargetEnd = getAdjacentEnd(seg2, _3PRIME);


	stPinchBlock *block = NULL;
	stSet *seen = stSet_construct();

	stList *stack = stList_construct();
	stPinchEnd *end = startEnd;
	stList_append(stack, end);
	while(stList_length(stack) > 0) {
		end = stList_pop(stack);
		block = stPinchEnd_getBlock(end);
		fprintf(stderr, "Arrived at %s\n", endInfo(end));

		stPinchEnd otherEnd = stPinchEnd_constructStatic(block, !stPinchEnd_getOrientation(end));

		if (leftTargetEnd && stPinchEnd_equalsFn(&otherEnd, leftTargetEnd)) {
			stList_destruct(stack);
			stPinchEnd *finalEnd = stPinchEnd_construct(stPinchSegment_getBlock(seg2), _3PRIME);
			return finalEnd;
		}
		else if (rightTargetEnd 
				&& stPinchEnd_equalsFn(&otherEnd, rightTargetEnd)) {
			stList_destruct(stack);
			stPinchEnd *finalEnd = stPinchEnd_construct(stPinchSegment_getBlock(seg2), _5PRIME);
			return finalEnd;
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
	return NULL;
}

bool pinchCreatesCycle(stPinchSegment *seg1, stPinchSegment *seg2, bool pinchOrientation) {
	stPinchEnd *path1End = directedWalk(seg1, seg2, _5PRIME);
	stPinchEnd *path2End = directedWalk(seg1, seg2, _3PRIME);

	//pinching these segments will create a cycle if there is a path 
	//connecting them, starting and ending at opposite ends of the newly
	//created block
	bool opposite = 
		(path1End && (pinchOrientation ^ (stPinchEnd_getOrientation(path1End) == _3PRIME))) ||
		(path2End && (pinchOrientation ^ (stPinchEnd_getOrientation(path2End) == _5PRIME)));
	return opposite;
}

bool acyclicFilterFn(stPinchSegment *seg1, stPinchSegment *seg2) {
	if (singleCopyFilterFn(seg1, seg2)) return true;
	
	return (pinchCreatesCycle(seg1, seg2, pinchOrientation));
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
		pinchOrientation = pinch->strand;

		stPinchThread_filterPinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand, acyclicFilterFn);

		if (!getOrdering(threadSet)) {
			fprintf(stderr, "Cycle created by previous pinch.\n");
			break;
		}


	}
	fclose(sequencesFile);
	stPinchIterator_destruct(pinchIterator);
	return threadSet;
}


