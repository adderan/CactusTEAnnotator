#include <stdio.h>
#include "bioioC.h"
#include "pairwiseAlignment.h"

#include "stPinchGraphs.h"
#include "stPinchIterator.h"

#include "repeatGraphs.h"


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

void printSegmentInfo(stPinchSegment *seg) {
	fprintf(stderr, "Arrived at segment %ld-%ld on thread %ld\n", 
				stPinchSegment_getStart(seg), 
				stPinchSegment_getStart(seg) + stPinchSegment_getLength(seg), 
				stPinchThread_getName(stPinchSegment_getThread(seg)));
}

stSet *getThreads(stPinchSegment *segment) {
	stSet *threads = stSet_construct();
	stPinchBlock *block = stPinchSegment_getBlock(segment);
	if (block) {
		stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
		stPinchSegment *otherSeg = NULL;
		while((otherSeg = stPinchBlockIt_getNext(&blockIt)) != NULL) {
			stSet_insert(threads, (void*) stPinchThread_getName(stPinchSegment_getThread(otherSeg)));
		}
	}
	else {
		stSet_insert(threads, (void*) stPinchThread_getName(stPinchSegment_getThread(segment)));
	}
	return threads;
}

bool singleCopyFilterFn(stPinchSegment *seg1, stPinchSegment *seg2) {
	if (stPinchSegment_getThread(seg1) == stPinchSegment_getThread(seg2)) return true;
	bool filter = false;

	stSet *threads1 = getThreads(seg1);
	stSet *threads2 = getThreads(seg2);

	stSet *intersect = stSet_getIntersection(threads1, threads2);
	if (stSet_size(intersect) > 0) filter = true;
	stSet_destruct(threads1);
	stSet_destruct(threads2);
	stSet_destruct(intersect);
	return filter;
}


void printBiedgedGraph(stPinchThreadSet *threadSet, stHash *coloring, FILE *gvizFile) {
	stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
	fprintf(gvizFile, "strict graph {\n");
	stPinchBlock *block;
	while((block = stPinchThreadSetBlockIt_getNext(&blockIt))) {
		stPinchEnd *end1 = stPinchEnd_construct(block, 0);
		stPinchEnd *end2 = stPinchEnd_construct(block, 1);

		//end nodes
		char *end1Color = (stHash_search(coloring, end1)) ? COLORNAME(stHash_search(coloring, end1)) : "grey";
		char *end2Color = (stHash_search(coloring, end2)) ? COLORNAME(stHash_search(coloring, end2)) : "grey";
		fprintf(gvizFile, "\t%ld[style=filled fillcolor=%s];\n", stPinchEnd_hashFn(end1), end1Color);
		fprintf(gvizFile, "\t%ld[style=filled fillcolor=%s];\n", stPinchEnd_hashFn(end2), end2Color);

				
		//block edge
		fprintf(gvizFile, "\t%ld -- %ld[dir=none color=blue penwidth=10];\n", stPinchEnd_hashFn(end1), stPinchEnd_hashFn(end2));

		//fprintf(gvizFile, "{rank = same; %ld; %ld;}\n", stPinchEnd_hashFn(end1), stPinchEnd_hashFn(end2));
		
		//near adjacencies
		stSet *nearAdjacent = stPinchEnd_getConnectedPinchEnds(end1);
		stSetIterator *nearAdjacentIt = stSet_getIterator(nearAdjacent);
		stPinchEnd *adjEnd;
		while((adjEnd = stSet_getNext(nearAdjacentIt))) {
			fprintf(gvizFile, "\t%ld -- %ld;\n", stPinchEnd_hashFn(end1), stPinchEnd_hashFn(adjEnd));

		}
		stSet_destructIterator(nearAdjacentIt);

		stSet *farAdjacent = stPinchEnd_getConnectedPinchEnds(end2);
		stSetIterator *farAdjacentIt = stSet_getIterator(farAdjacent);
		while((adjEnd = stSet_getNext(farAdjacentIt))) {
			fprintf(gvizFile, "\t%ld -- %ld;\n", stPinchEnd_hashFn(end2), stPinchEnd_hashFn(adjEnd));
		}
	}
	fprintf(gvizFile, "}\n");
}

/*If the graph is not acyclic, returns NULL. If the graph is acyclic, 
* return an ordering of the nodes such that all nodes reachable 
* from any given node occur after it in the ordering.*/
stList *getComponentOrdering(stPinchBlock *startBlock, stHash *coloring) {
	stList *stack = stList_construct();
	stSet *seen = stSet_construct();
	
	stList *ordering = stList_construct();

	stPinchEnd *leftEnd = stPinchEnd_construct(startBlock, _5PRIME);

	stList_append(stack, leftEnd);
	stHash_insert(coloring, leftEnd, BLACK);

	stPinchBlock *block;
	stPinchEnd *end;
	while (stList_length(stack) > 0) {
		end = stList_pop(stack);
		//fprintf(stderr, "Arrived at %s\n", endInfo(end));
		block = stPinchEnd_getBlock(end);

		stPinchEnd *oppositeEnd = stPinchEnd_construct(block, !stPinchEnd_getOrientation(end));

		assert(stHash_search(coloring, end));
		if (stHash_search(coloring, end) ==
				stHash_search(coloring, oppositeEnd)) {
			//encountered cycle
			fprintf(stderr, "Encountered cycle at end %ld, both ends %p.\n", stPinchEnd_hashFn(end), stHash_search(coloring, end));
			stList_destruct(stack);
			stList_destruct(ordering);
			return NULL;
		}
		stHash_insert(coloring, oppositeEnd, OP(stHash_search(coloring, end)));

		if (!stSet_search(seen, block)) {
			stSet_insert(seen, block);
			stSet *nearAdjacentEnds = stPinchEnd_getConnectedPinchEnds(end);
			//fprintf(stderr, "Found %ld connected ends on near side\n", stSet_size(nearAdjacentEnds));
			stSetIterator *adjEndIt = stSet_getIterator(nearAdjacentEnds);
			stPinchEnd *adjEnd;
			while((adjEnd = stSet_getNext(adjEndIt))) {
				void *colorToTagEnd = OP(stHash_search(coloring, end));
				if (stHash_search(coloring, adjEnd) == OP(colorToTagEnd)) {
					fprintf(stderr, "About to encounter cycle\n");
				}
				stHash_insert(coloring, adjEnd, colorToTagEnd);
				stList_append(stack, adjEnd);
			}

			stSet_destructIterator(adjEndIt);
			stSet *farAdjacentEnds = stPinchEnd_getConnectedPinchEnds(oppositeEnd);
			//fprintf(stderr, "Found %ld connected ends on far side.\n", stSet_size(farAdjacentEnds));
			adjEndIt = stSet_getIterator(farAdjacentEnds);
			while((adjEnd = stSet_getNext(adjEndIt))) {
				void *colorToTagEnd = OP(stHash_search(coloring, oppositeEnd));
				if (stHash_search(coloring, adjEnd) == OP(colorToTagEnd)) {
					fprintf(stderr, "About to encounter cycle\n");
				}
				stHash_insert(coloring, adjEnd, colorToTagEnd);
				stList_append(stack, adjEnd);
			}
			stSet_destructIterator(adjEndIt);
		}
	}
	stList_destruct(stack);
	return ordering;
}

stList *getOrdering(stPinchThreadSet *graph, stHash *coloring) {
	stSortedSet *threadComponents = stPinchThreadSet_getThreadComponents(graph);
	fprintf(stderr, "Found %ld components\n", stSortedSet_size(threadComponents));
	stList *ordering = stList_construct();
	if (!coloring) {
		coloring = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, (void(*)(void*)) stPinchEnd_destruct, NULL);
	}

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
		stList *componentOrdering = getComponentOrdering(block, coloring);
		if (!componentOrdering) {
			ordering = NULL;
			goto out;
		}
		stList_appendAll(ordering, componentOrdering);

	}
out:
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

//inefficient function for checking whether there is a directed
//walk (a path of alternating block and adjacency edges) between
//two segments.
bool directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool startDirection) {
	assert(stPinchSegment_getThread(seg1) != stPinchSegment_getThread(seg2));

	stPinchSegment *seg = seg1;
	stPinchBlock *block = stPinchSegment_getBlock(seg);


	stPinchEnd *leftTargetEnd = (block) ? stPinchEnd_construct(stPinchSegment_getBlock(seg2), _5PRIME) : NULL;
	stPinchEnd *rightTargetEnd = (block) ? stPinchEnd_construct(stPinchSegment_getBlock(seg2), _3PRIME) : NULL;

	stSet *seen = stSet_construct();

	stList *stack = stList_construct();
	stList *directions = stList_construct();
	stList_append(stack, seg);
	stList_append(directions, (void*)startDirection);
	bool curDirection = startDirection;
	bool reachable = false;
	while(seg || stList_length(stack) > 0) {
		if (!seg) {
			seg = stList_pop(stack);
			curDirection = stList_pop(directions);
		}
		if (seg == seg2) {
			reachable = true;
			goto out;
		}

		block = stPinchSegment_getBlock(seg);

		if (block && leftTargetEnd && rightTargetEnd) {
			stPinchEnd *leftEnd = stPinchEnd_construct(block, 0);
			stPinchEnd *rightEnd = stPinchEnd_construct(block, 1);
			if (stPinchEnd_equalsFn(leftEnd, leftTargetEnd) || 
					stPinchEnd_equalsFn(rightEnd, rightTargetEnd) || 
					stPinchEnd_equalsFn(rightEnd, leftTargetEnd) || 
					stPinchEnd_equalsFn(leftEnd, rightTargetEnd)) {
				reachable = true;
				goto out;
			}
			stPinchEnd_destruct(leftEnd);
			stPinchEnd_destruct(rightEnd);
		}
		if (block && !stSet_search(seen, block)) {
			stSet_insert(seen, block);
			stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
			stPinchSegment *nextSeg;
			while((nextSeg = stPinchBlockIt_getNext(&blockIt))) {
				if (nextSeg != seg) {
					stList_append(stack, nextSeg);

					//switch directions if the next segment 
					//or the current segment is reversed relative 
					//to the block
					bool nextSegDirection = curDirection 
						^ (stPinchSegment_getBlockOrientation(nextSeg)==0)
						^ (stPinchSegment_getBlockOrientation(seg)==0);
					stList_append(directions, (void*)nextSegDirection);
				}
			}

		}
		seg = (curDirection == _5PRIME) ? 
			stPinchSegment_get5Prime(seg) : stPinchSegment_get3Prime(seg);

	}
out:
	stList_destruct(stack);
	stList_destruct(directions);
	stSet_destruct(seen);
	return reachable;
}


bool acyclicFilterFn(stPinchSegment *seg1, stPinchSegment *seg2) {
	if (singleCopyFilterFn(seg1, seg2)) return true;
	
	return (directedWalk(seg1, seg2, _5PRIME) || directedWalk(seg1, seg2, _3PRIME));
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

		stPinchThread_filterPinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand, acyclicFilterFn);

		fprintf(stderr, "Checking pinch\n");
		if (!getOrdering(threadSet, NULL)) {
			fprintf(stderr, "Cycle created by previous pinch.\n");
			stHash *coloring = stHash_construct3(stPinchEnd_hashFn, 
					stPinchEnd_equalsFn, 
					(void(*)(void*)) stPinchEnd_destruct, NULL);

			getOrdering(threadSet, coloring);
			FILE *afterPinchFile = fopen("after_bad_pinch.gvz", "w");
			printBiedgedGraph(threadSet, coloring, afterPinchFile);
			fclose(afterPinchFile);
			stHash_destruct(coloring);

			stPinchUndo *undo = stPinchThread_prepareUndo(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand);
			assert(undo);
			stPinchThreadSet_undoPinch(threadSet, undo);
			stHash *coloring2 = stHash_construct3(stPinchEnd_hashFn, 
					stPinchEnd_equalsFn, 
					(void(*)(void*)) stPinchEnd_destruct, NULL);
			getOrdering(threadSet, coloring2);
			FILE *beforePinchFile = fopen("before_bad_pinch.gvz", "w");
			printBiedgedGraph(threadSet, coloring2, beforePinchFile);
			fclose(beforePinchFile);
			stHash_destruct(coloring2);
			break;
		}
	}
	fclose(sequencesFile);
	stPinchIterator_destruct(pinchIterator);
	return threadSet;
}
