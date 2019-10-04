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

stSortedSet *getThreads(stPinchSegment *segment) {
	stSortedSet *threads = stSortedSet_construct();
	stPinchBlock *block = stPinchSegment_getBlock(segment);
	if (block) {
		stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
		stPinchSegment *otherSeg = NULL;
		while((otherSeg = stPinchBlockIt_getNext(&blockIt)) != NULL) {
			stSortedSet_insert(threads, (void*) stPinchThread_getName(stPinchSegment_getThread(otherSeg)));
		}
	}
	else {
		stSortedSet_insert(threads, (void*) stPinchThread_getName(stPinchSegment_getThread(segment)));
	}
	return threads;
}

bool singleCopyFilterFn(stPinchSegment *seg1, stPinchSegment *seg2) {
	if (stPinchSegment_getThread(seg1) == stPinchSegment_getThread(seg2)) return true;
	bool filter = false;

	stSortedSet *threads1 = getThreads(seg1);
	stSortedSet *threads2 = getThreads(seg2);

	stSortedSet *intersect = stSortedSet_getIntersection(threads1, threads2);
	if (stSortedSet_size(intersect) > 0) filter = true;
	stSortedSet_destruct(threads1);
	stSortedSet_destruct(threads2);
	stSortedSet_destruct(intersect);
	return filter;
}

void printBiedgedGraph(stPinchThreadSet *threadSet, char *gvizFilename) {
	stHash *coloring = stHash_construct3(stPinchEnd_hashFn, 
				stPinchEnd_equalsFn, 
				(void(*)(void*)) stPinchEnd_destruct, NULL);
	getOrdering(threadSet, coloring);

	FILE *gvizFile = fopen(gvizFilename, "w");
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
	stHash_destruct(coloring);
	fprintf(gvizFile, "}\n");
	fclose(gvizFile);
}

stList *getOrdering(stPinchThreadSet *threadSet, stHash *coloring) {
	stList *ordering = stList_construct();
	stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIterator(threadSet);
	//if all the blocks with incoming adjacencies have been seen, this 
	//block can be added to the ordering
	stSet *incomingEnds = stPinchEnd_getConnectedPinchEnds(end);
	stSetIterator *incomingEndIt = stSet_getIterator(incomingEnds);
	stPinchEnd *incomingEnd;
	bool addToOrdering = true;
	while((incomingEnd = stSet_getNext(incomingEndIt)) != NULL) {
		if (!stSet_search(seen, stPinchEnd_getBlock(incomingEnd))) {
			addToOrdering = false;
			break;
		}
	}

}

/*If the graph is not acyclic, returns NULL. If the graph is acyclic, 
* return a dictionary assigning an orientation to each end.*/
bool getColoring2(stPinchBlock *startBlock, stHash *coloring) {
	stList *stack = stList_construct();
	stSet *seen = stSet_construct();

	stPinchEnd *leftEnd = stPinchEnd_construct(startBlock, _5PRIME);

	stList_append(stack, leftEnd);
	stHash_insert(coloring, leftEnd, BLACK);

	stPinchBlock *block;
	stPinchEnd *end;
	while (stList_length(stack) > 0) {
		end = stList_pop(stack);
		block = stPinchEnd_getBlock(end);

		stPinchEnd *oppositeEnd = stPinchEnd_construct(block, !stPinchEnd_getOrientation(end));

		assert(stHash_search(coloring, end));
		if (stHash_search(coloring, end) ==
				stHash_search(coloring, oppositeEnd)) {
			//encountered cycle
			fprintf(stderr, "Encountered cycle at end %ld, both ends %p.\n", stPinchEnd_hashFn(end), stHash_search(coloring, end));
			stList_destruct(stack);
			return false;
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
					return NULL;
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
					return NULL;
				}
				stHash_insert(coloring, adjEnd, colorToTagEnd);
				stList_append(stack, adjEnd);
			}
			stSet_destructIterator(adjEndIt);
		}
	}
	stList_destruct(stack);
	return true;
}

stList *getColoring(stPinchThreadSet *graph) {
	stHash *coloring = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, (void(*)(void*)) stPinchEnd_destruct, NULL);

	stSortedSet *threadComponents = stPinchThreadSet_getThreadComponents(graph);

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
		bool acyclic = getColoring2(block, coloring);
		if (!acyclic) {
			stHash_destruct(coloring);
			coloring = NULL;
			goto out;
		}

	}
out:
	stSortedSet_destructIterator(componentsIt);
	return coloring;
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
	stPinchBlock *block;
	stSet *seenBlocks = stSet_construct();
	stList *stack = stList_construct();
	stList *directions = stList_construct();
	stPinchSegment *segment = seg1;
	stList_append(stack, segment);
	stList_append(directions, (void*)startDirection);
	bool curDirection = startDirection;
	bool walk = false;
	while(segment || stList_length(stack) > 0) {
		if (!segment) {
			segment = stList_pop(stack);
			curDirection = stList_pop(directions);
		}
		if (segment == seg2) {
			walk = true;
			goto out;
		}
		block = stPinchSegment_getBlock(segment);
		if (block && !stSet_search(seenBlocks, block)) {
			stSet_insert(seenBlocks, block);
			stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
			stPinchSegment *nextSeg;
			while((nextSeg = stPinchBlockIt_getNext(&blockIt))) {
				if (nextSeg != segment) {
					stList_append(stack, nextSeg);

					//switch directions if the next segment 
					//or the current segment is reversed relative 
					//to the block
					bool nextSegDirection = curDirection 
						^ (stPinchSegment_getBlockOrientation(nextSeg)==0)
						^ (stPinchSegment_getBlockOrientation(segment)==0);
					stList_append(directions, (void*)nextSegDirection);
				}
			}

		}
		segment = (curDirection == _5PRIME) ? 
			stPinchSegment_get5Prime(segment) : stPinchSegment_get3Prime(segment);

	}
out:
	stList_destruct(stack);
	stList_destruct(directions);
	stSet_destruct(seenBlocks);
	return walk;
}

bool acyclicFilterFn(stPinchSegment *seg1, stPinchSegment *seg2) {
	if (singleCopyFilterFn(seg1, seg2)) return true;
	return (directedWalk(seg1, seg2, _3PRIME) || directedWalk(seg1, seg2, _5PRIME));
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

stPinchThreadSet *initializeGraph(char *sequencesFilename) {
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
	return threadSet;
}

void applyPinch(stPinchThreadSet *threadSet, stPinch *pinch) {
	if (pinch->strand == 0) return;
	stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
	stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
	stPinchThread_filterPinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand, acyclicFilterFn);
}


stPinchThreadSet *buildRepeatGraph(char *sequencesFilename, char *alignmentsFilename) {

	stPinchThreadSet *threadSet = initializeGraph(sequencesFilename);
#ifndef NDEBUG
	//maintain another copy of the graph to keep track of what it looked
	//like before applying a bad pinch
	stPinchThreadSet *threadSet_debug = initializeGraph(sequencesFilename);
#endif
	
	stPinchIterator *pinchIterator = stPinchIterator_constructFromFile(alignmentsFilename);
	stPinch *pinch = NULL;

	while((pinch = stPinchIterator_getNext(pinchIterator)) != NULL) {
		applyPinch(threadSet, pinch);

#ifndef NDEBUG
		stList *ordering = getOrdering(threadSet, NULL);
		if (!ordering) {
			fprintf(stderr, "Pinch created cycle.\n");
			printBiedgedGraph(threadSet, "after_bad_pinch.gvz");
			printBiedgedGraph(threadSet_debug, "before_bad_pinch.gvz");
			exit(1);
		}
		applyPinch(threadSet_debug, pinch);
#endif
	}
	assert(getOrdering(threadSet, NULL));
	stPinchIterator_destruct(pinchIterator);
	return threadSet;
}

//Finds the highest-weight traversal through a connected thread component
//in an acyclic pinch graph, given an ordered listing of the nodes in 
//the component
char *getConsensusPath(stPinchThreadSet *graph, stHash *coloring) {
	char *path = malloc(sizeof(char)*stList_length(component));
	fprintf(stderr, "Component with %ld blocks\n", stList_length(component));
	stHash *blockIndex = stHash_construct();

	for (int i = 0; i < stList_length(component); i++) {
		stPinchBlock *block = stList_get(component, i);
		stPinchEnd *rightEnd = stPinchEnd_construct(block, _3PRIME);
		stPinchEnd *leftEnd = stPinchEnd_construct(block, _5PRIME);
		stPinchEnd *blackEnd = (stHash_search(coloring, leftEnd) == BLACK) ? leftEnd : rightEnd;


		stList *precedingEnds = stPinchEnd_getConnectedPinchEnds(blackEnd);
		for (int j = 0; j < stList_length(precedingEnds); j++) {
			stPinchBlock *precedingBlock = stPinchEnd_getBlock(stList_get(precedingEnds, j));
			assert(stHash_search(blockIndex, precedingBlock) != NULL);
			assert(stHash_search(blockIndex, precedingBlock) < i);
		}
	}

	return NULL;
}
