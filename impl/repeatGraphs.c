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

void printBiedgedGraph(stPinchThreadSet *graph, char *gvizFilename) {
	stHash *coloring = graphIsAcyclic(graph);

	FILE *gvizFile = fopen(gvizFilename, "w");
	stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(graph);
	fprintf(gvizFile, "strict graph {\n");
	stPinchBlock *block;
	while((block = stPinchThreadSetBlockIt_getNext(&blockIt))) {
		stPinchEnd *end1 = stPinchEnd_construct(block, _3PRIME);
		stPinchEnd *end2 = stPinchEnd_construct(block, _5PRIME);

		stPinchEnd *blackEnd = (stHash_search(coloring, end1) == BLACK) ? end1 : end2;
		stPinchEnd *redEnd = (blackEnd == end1) ? end2 : end1;

		//end nodes
		fprintf(gvizFile, "\t%ld[style=filled fillcolor=grey];\n", stPinchEnd_hashFn(blackEnd));
		fprintf(gvizFile, "\t%ld[style=filled fillcolor=red];\n", stPinchEnd_hashFn(redEnd));

				
		//block edge
		fprintf(gvizFile, "\t%ld -- %ld[dir=none color=blue penwidth=10];\n", stPinchEnd_hashFn(blackEnd), stPinchEnd_hashFn(redEnd));

		//fprintf(gvizFile, "{rank = same; %ld; %ld;}\n", stPinchEnd_hashFn(end1), stPinchEnd_hashFn(end2));
		
		//near adjacencies
		stSet *nearAdjacent = stPinchEnd_getConnectedPinchEnds(blackEnd);
		stSetIterator *nearAdjacentIt = stSet_getIterator(nearAdjacent);
		stPinchEnd *adjEnd;
		while((adjEnd = stSet_getNext(nearAdjacentIt))) {
			fprintf(gvizFile, "\t%ld -- %ld;\n", stPinchEnd_hashFn(blackEnd), stPinchEnd_hashFn(adjEnd));

		}
		stSet_destructIterator(nearAdjacentIt);

		stSet *farAdjacent = stPinchEnd_getConnectedPinchEnds(redEnd);
		stSetIterator *farAdjacentIt = stSet_getIterator(farAdjacent);
		while((adjEnd = stSet_getNext(farAdjacentIt))) {
			fprintf(gvizFile, "\t%ld -- %ld;\n", stPinchEnd_hashFn(redEnd), stPinchEnd_hashFn(adjEnd));
		}
	}
	stHash_destruct(coloring);
	fprintf(gvizFile, "}\n");
	fclose(gvizFile);
}

stList *getOrdering(stPinchBlock *startBlock) {
	stList *ordering = stList_construct();
	stList *stack = stList_construct();

	stSet *seen = stSet_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, (void*) stPinchEnd_destruct);

	stList_append(stack, startBlock);
	stPinchBlock *block;
	stPinchEnd *end;
	while(stList_length(stack) > 0) {
		end = stList_pop(stack);
		block = stPinchEnd_getBlock(end);


		stSet *incomingAdjacencies = stPinchEnd_getConnectedPinchEnds(end);
		stSetIterator *incomingEndIt = stSet_getIterator(incomingAdjacencies);
		stPinchEnd *incomingEnd;
		bool predecessorsSeen = true;
		while((incomingEnd = stSet_getNext(incomingEndIt))) {
			if (!stSet_search(seen, incomingEnd)) {
				predecessorsSeen = false;
				//jump to the other end of the block, since we're going
				//backwards now
				stPinchBlock *incomingBlock = stPinchEnd_getBlock(incomingEnd);
				stPinchEnd *blackEnd = stPinchEnd_construct(incomingBlock, !stPinchEnd_getOrientation(incomingEnd));
				stList_append(stack, blackEnd);
			}
		}
		if (predecessorsSeen) {
			stList_append(ordering, end);
			stPinchEnd *oppositeEnd = stPinchEnd_construct(block, !stPinchEnd_getOrientation(end));
			stSet *outgoingAdjacencies = stPinchEnd_getConnectedPinchEnds(oppositeEnd);
			stSetIterator *outgoingEndIt = stSet_getIterator(outgoingAdjacencies);
			stPinchEnd *outgoingEnd;
			while((outgoingEnd = stSet_getNext(outgoingEndIt))) {
				stList_append(stack, outgoingEnd);
			}
			stPinchEnd_destruct(oppositeEnd);
		}

	}
	return ordering;
}

bool graphIsAcyclic(stPinchThreadSet *graph) {
	stList *threadComponents = stPinchThreadSet_getThreadComponents(graph);

	stListIterator *componentsIt = stList_getIterator(threadComponents);
	stList *component;
	bool acyclic = true;
	while((component = stList_getNext(componentsIt))) {
		stPinchThread *thread = stList_peek(component);
		stPinchBlock *startBlock = getFirstBlock(thread);
		if (!componentIsAcyclic(startBlock)) {
			acyclic = false;
			goto out;
		}
	}
out:
	stList_destruct


}


/*If the graph is not acyclic, returns NULL. If the graph is acyclic, 
* return a list of edges in the DAG corresponding to this pinch graph,
* in ordering implied by traversing the graph. */
bool componentIsAcyclic(stPinchBlock *startBlock) {
	stList *stack = stList_construct();

	stSet *seen = stSet_construct();

	stHash *coloring = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, (void(*)(void*)) stPinchEnd_destruct, NULL);

	stPinchEnd *leftEnd = stPinchEnd_construct(startBlock, _5PRIME);

	stList_append(stack, leftEnd);
	stHash_insert(coloring, leftEnd, BLACK);
	
	stPinchBlock *block;
	stPinchEnd *end;
	while (stList_length(stack) > 0) {
		end = stList_pop(stack);
		block = stPinchEnd_getBlock(end);
		fprintf(stderr, "Arrived at block %p\n", block);

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
			//fprintf(stderr, "Found %ld connected ends on near side\n", stSet_size(nearAdjacentEnds));
			stSet *incomingEnds = stPinchEnd_getConnectedPinchEnds(end);
			stSetIterator *incomingEndIt = stSet_getIterator(incomingEnds);
			stPinchEnd *adjEnd;
			while((adjEnd = stSet_getNext(incomingEndIt))) {
				void *colorToTagEnd = OP(stHash_search(coloring, end));
				if (stHash_search(coloring, adjEnd) == OP(colorToTagEnd)) {
					fprintf(stderr, "About to encounter cycle\n");
					return NULL;
				}
				stHash_insert(coloring, adjEnd, colorToTagEnd);
				stList_append(stack, adjEnd);
			}

			stSet_destructIterator(incomingEndIt);
			stSet *outgoingEnds = stPinchEnd_getConnectedPinchEnds(oppositeEnd);
			stSetIterator *outgoingEndIt = stSet_getIterator(outgoingEnds);
			outgoingEndIt = stSet_getIterator(outgoingEnds);
			while((adjEnd = stSet_getNext(outgoingEndIt))) {
				void *colorToTagEnd = OP(stHash_search(coloring, oppositeEnd));
				if (stHash_search(coloring, adjEnd) == OP(colorToTagEnd)) {
					fprintf(stderr, "About to encounter cycle\n");
					return NULL;
				}
				stHash_insert(coloring, adjEnd, colorToTagEnd);
				stList_append(stack, adjEnd);
			}
			stSet_destructIterator(outgoingEndIt);
		}
	}
	stList_destruct(stack);
	stSet_destruct(seen);
	stSet_destruct(seenEnds);
	return coloring;
}

stPinchBlock *getFirstBlock(stPinchThread *thread) {
	stPinchSegment *segment = stPinchThread_getFirst(thread);
	while(segment && (stPinchSegment_getBlock(segment) == NULL)) {
		segment = stPinchSegment_get3Prime(segment);
	}
	if (!segment) {
		segment = stPinchThread_getFirst(thread);
		while(segment && (stPinchSegment_getBlock(segment) == NULL)) {
			segment = stPinchSegment_get3Prime(segment);
		}
	}
	if (!segment) return NULL;
	return stPinchSegment_getBlock(segment);
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

stPinchThreadSet *initializeGraph(stHash *sequences) {
	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	stHashIterator *sequencesIt = stHash_getIterator(sequences);
	void *threadName;
	while((threadName = stHash_getNext(sequencesIt))) {
		char *sequence = stHash_search(sequences, threadName);
		stPinchThreadSet_addThread(threadSet, (int64_t) threadName, 0, strlen(sequence));
	}
	return threadSet;
}

void applyPinch(stPinchThreadSet *threadSet, stPinch *pinch) {
	if (pinch->strand == 0) return;
	stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
	stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
	stPinchThread_filterPinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand, acyclicFilterFn);
}

stPinchThreadSet *buildRepeatGraph(stHash *sequences, char *alignmentsFilename) {

	stPinchThreadSet *threadSet = initializeGraph(sequences);
#ifndef NDEBUG
	//maintain another copy of the graph to keep track of what it looked
	//like before applying a bad pinch
	stPinchThreadSet *threadSet_debug = initializeGraph(sequences);
#endif
	
	stPinchIterator *pinchIterator = stPinchIterator_constructFromFile(alignmentsFilename);
	stPinch *pinch = NULL;

	while((pinch = stPinchIterator_getNext(pinchIterator)) != NULL) {
		applyPinch(threadSet, pinch);

#ifndef NDEBUG
		/*
		stList *ordering = getOrdering(threadSet);
		if (!ordering) {
			fprintf(stderr, "Pinch created cycle.\n");
			printBiedgedGraph(threadSet, "after_bad_pinch.gvz");
			printBiedgedGraph(threadSet_debug, "before_bad_pinch.gvz");
			exit(1);
		}
		*/
		applyPinch(threadSet_debug, pinch);
#endif
	}
	stPinchIterator_destruct(pinchIterator);
	return threadSet;
}

stPinchSegment *getSegment(stPinchBlock *block, int64_t threadName) {
	stPinchSegment *segment;
	stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
	while((segment = stPinchBlockIt_getNext(&blockIt))) {
		if (stPinchThread_getName(stPinchSegment_getThread(segment)) == threadName) break;
	}
	return segment;
}

stList *traversePath(stPinchThreadSet *threadSet, stList *endsInPath, stHash *sequences) {
	stList *path = stList_construct();

	for (int64_t i = 0; i < stList_length(endsInPath); i++) {
		if (i > 0) {
			//Fill in the adjacency connecting the previous block
			//to this one
			stPinchEnd *prevEnd = stList_get(endsInPath, i - 1);
			stPinchEnd *end1 = stPinchEnd_construct(stPinchEnd_getBlock(prevEnd), !stPinchEnd_getOrientation(prevEnd));
			stPinchEnd *end2 = stList_get(endsInPath, i);
			assert(end1);
			assert(end2);

			stPinchBlock *block1 = stPinchEnd_getBlock(end1);
			stPinchBlock *block2 = stPinchEnd_getBlock(end2);
			stSortedSet *threads1 = 
				getThreads(stPinchBlock_getFirst(block1));

			stSortedSet *threads2 = 
				getThreads(stPinchBlock_getFirst(block2));

			stSortedSet *sharedThreadsSet = stSortedSet_getIntersection(threads1, threads2);
			stList *sharedThreads = stSortedSet_getList(sharedThreadsSet);

			stList *candidateAdjacencySeqs = stList_construct();
			for (int j = 0; j < stList_length(sharedThreads); j++) {
				int64_t threadName = (int64_t) stList_get(sharedThreads, j);

				stPinchSegment *seg1 = getSegment(block1, threadName);
				stPinchSegment *seg2 = getSegment(block2, threadName);
				
				int64_t adjacencyLength = stPinchSegment_getStart(seg2) - stPinchSegment_getStart(seg1) - stPinchSegment_getLength(seg1);
				assert(adjacencyLength > 0);

				char *adj = malloc(sizeof(char)* (adjacencyLength + 1));
				strncpy(stHash_search(sequences, (void*)threadName), adj, adjacencyLength);
				adj[adjacencyLength] = '\0';
				stList_append(candidateAdjacencySeqs, adj);
				
			}

			//pick the shortest one
			int64_t minLength = INT_MAX;
			char *adjacencySequence = NULL;
			for (int j = 0; j < stList_length(candidateAdjacencySeqs); j++) {
				char *adj = stList_get(candidateAdjacencySeqs, j);
				if (strlen(adj) < minLength) {
					minLength = strlen(adj);
					free(adjacencySequence);
					adjacencySequence = adj;
				}
				else {
					free(adj);
				}
			}
			assert(adjacencySequence);

			stList_append(path, adjacencySequence);

		}

		stPinchEnd *end = stList_get(endsInPath, i);
		assert(end);
		stPinchBlock *block = stPinchEnd_getBlock(end);
		stPinchSegment *segment = stPinchBlock_getFirst(block);
		int64_t threadName = stPinchThread_getName(stPinchSegment_getThread(segment));
		char *blockSequence = malloc(sizeof(char) * (stPinchSegment_getLength(segment) + 1));
		strncpy(stHash_search(sequences, (void*)threadName), blockSequence, stPinchSegment_getLength(segment));
		blockSequence[stPinchSegment_getLength(segment) + 1] = '\0';

		stList_append(path, blockSequence);
	}
	return path;
}

stList *heaviestPath(stPinchThreadSet *graph, stList *ordering) {
	int N = stList_length(ordering);
	int64_t *paths = calloc(sizeof(int), N);
	int64_t *score = calloc(sizeof(int), N);

	stHash *blockIndex = stHash_construct();
	for (int64_t i = 0; i < N; i++) {
		stHash_insert(blockIndex, stList_get(ordering, i), (void*)i);
	}

	int64_t bestScore = 0;
	int64_t bestEndpoint = -1;

	for (int64_t i = N - 1; i >= 0; i--) {
		stPinchEnd *end = stList_get(ordering, i);
		stPinchBlock *block = stPinchEnd_getBlock(end);

		int bestRight = -1;
		int maxWeight = 0;

		stSet *precedingEnds = stPinchEnd_getConnectedPinchEnds(end);
		stSetIterator *endIt = stSet_getIterator(precedingEnds);
		stPinchEnd *precedingEnd;
		while((precedingEnd = stSet_getNext(endIt))) {
			stPinchBlock *precedingBlock = stPinchEnd_getBlock(precedingEnd);
			int64_t weight = stPinchBlock_getDegree(precedingBlock) * stPinchBlock_getLength(precedingBlock);
			if (weight > maxWeight) {
				bestRight = (int64_t) stHash_search(blockIndex, block);
				maxWeight = weight;
			}

		}

		paths[i] = bestRight;

		if (score[i] > bestScore) {
			bestEndpoint = i;
			bestScore = score[i];
		}
	}
	
	//traceback
	stList *path = stList_construct();
	for (int i = bestEndpoint; i >= 0; i = paths[i]) {
		assert(stList_get(ordering, i));
		stList_append(path, stList_get(ordering, i));
	}

	return path;
}
